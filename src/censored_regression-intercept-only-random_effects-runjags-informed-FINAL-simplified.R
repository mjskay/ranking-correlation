library(coda)
library(rjags)
library(runjags)
library(survival)
library(gamlss)
library(gamlss.cens)
library(ggplot2)
library(grid)
library(magrittr)
library(plyr)
library(dplyr)
library(tidyr)
library(tidybayes)
library(stringr)
library(scales)

memory.limit(8000)

source("src/openGraphSaveGraph.R")

source("src/clean-data.R")

#------------------------------------------------------------------------------
# THE MODEL.
include_predictions = TRUE
final_model = FALSE
include_typical = TRUE

modelstring = paste("
model {
	#MODEL
	#core model
	for (i in 1:n) {
		# latent variable log-linear model
		mu[i] <- b[visandsign[i],1] + b[visandsign[i],2]*r[i] + 
		         b[visandsign[i],3]*approach_value[i] + b[visandsign[i],4]*approach_value[i]*r[i] +
				 u[participant[i]]
    	y[i] ~ dlnorm(mu[i], tau[visandsign[i]])

		# right-censoring at this value of r
		censored[i] ~ dinterval(y[i], censoring_threshold[i])
  	}
    
	#participant random effects
	for (p in 1:n_participant) {
		u[p] ~ dnorm(0, u_tau[participant_visandsign[p]])
	}

	#priors
	for (v in 1:n_visandsign) {
		#priors for coefficents of fixed effects
		b[v,1] ~ dnorm(log(.45), 1)	#prior on intercept is chance
		b[v,2] ~ dnorm(0, 0.05)
		b[v,3] ~ dnorm(0, 4)
		b[v,4] ~ dnorm(0, 4)

    	#prior on variance of latent variable
    	tau[v] ~ dgamma(1, 1)

		#prior on random effects
		u_tau[v] ~ dgamma(1, 1)
	}

	#EXPECTED PERFORMANCE ON RANDOM DATASET
	", if (include_typical) "
	typical_r ~ dunif(0.3,0.8)
	for (v in 1:n_visandsign) {
		typical_mu[v] <- b[v,1] + b[v,2]*typical_r
	}" else "", "

	#PREDICTIONS
	", if (include_predictions) "
	for (i in 1:n_pred) {
		pred_u[i] ~ dnorm(0, u_tau[pred_visandsign[i]])
		pred_mu[i] <- b[pred_visandsign[i],1] + b[pred_visandsign[i],2]*pred_r[i] +
					pred_u[i]
    	pred_y[i] ~ dlnorm(pred_mu[i], tau[pred_visandsign[i]])
	}" else "", "

}")
writeLines(modelstring,con="model.txt")


#------------------------------------------------------------------------------
# THE DATA.

#determine the visandsign each participant was assigned to
#N.B. this relies on the assumption (correct for the Harrison data)
#that all participants were assigned to only one visandsign
participant_visandsign = df %>%
    group_by(participant) %>%
    arrange(participant) %>%		#put in order so that indices line up in JAGS
    slice(1) %$%
    as.numeric(visandsign)

#build data list
data_list = df %>%
    mutate(y = ifelse(censored, NA, jnd)) %>%	#values of latent variable are observed data only if not censored
    select(y, r, approach_value, visandsign, participant, censoring_threshold, censored) %>%
    compose_data(participant_visandsign)

#------------------------------------------------------------------------------
# DATA FOR PREDICTIONS
if (include_predictions) {
    #add predictions to data list
    pred_df = expand.grid(
        pred_visandsign=1:length(levels(df$visandsign)), 
        pred_r=unique(df$r))
    data_list = data_list %>% 
        compose_data(pred=pred_df)
    #make pred_df more useful for later
    pred_df = pred_df %>%
        mutate(i=1:nrow(.), visandsign=factor(pred_visandsign, labels=levels(df$visandsign))) %>%    
        select(r=pred_r, visandsign, i)
}

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
#get maximum likelihood estimates of model parameters
coefs = ddply(df, ~ visandsign, function(df) {
        m = gamlss(Surv(censored_jnd, !censored) ~ r * approach_value, family=cens(LOGNO), data=df)
        c(
            sigma = exp(coef(m, "sigma")), 
            coef(m)
        )
    })

#y initialization should provide values for the latent variable only for observations 
#above the censoring threshold (values below the threshold are observed data so
#cannot be initialized here)
y_init = with(df, ifelse(censored, jnd, NA))

inits_list = list(
	b=as.matrix(coefs[,-1:-2]),	#first two columns are visandsign and tau
    tau=1/(coefs$sigma^2),
    y=y_init
)

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c("b" , "tau", "u_tau")  # The parameter(s) to be monitored.
if (include_predictions) parameters = c(parameters, "pred_y") 
if (include_typical) parameters = c(parameters, "typical_r", "typical_mu") 

if (!final_model) {
    jagsModel = run.jags("model.txt", data=data_list, monitor=parameters, initlist=inits_list, 
        method="parallel")
} else {
#    jagsModel = autorun.jags("model.txt", data=data_list, monitor=parameters, initlist=inits_list,
#        method="parallel", thin.sample=TRUE)    
    jagsModel = run.jags("model.txt", data=data_list, monitor=parameters, initlist=inits_list,
	    adapt=5000, burnin = 100000, sample = 10000, thin=10,
        method="parallel")
}


#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

codaSamples = as.mcmc.list(jagsModel)

checkConvergence = FALSE
if ( checkConvergence ) {
    plot(jagsModel, vars="^b", file="output/model-params-b.pdf")
    plot(jagsModel, vars="^tau", file="output/model-params-tau.pdf")
    plot(jagsModel, vars="^u_tau", file="output/model-params-u_tau.pdf")
    plot(jagsModel, vars="^typ", file="output/model-params-typical_mu.pdf")
    plot(jagsModel, vars="^pred", file="output/model-params-pred.pdf")

    summary(jagsModel)
    pdf(file="output/model-autocorr.pdf")
    autocorr.plot(codaSamples, ask=FALSE)
    dev.off()
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )
rm("codaSamples")
rm("jagsModel")

# Extract chain values:
bdf = extract_samples(mcmcChain, ~ b[visandsign, ..], df)
bdf_m = bdf %>%
    group_by(visandsign) %>%
    summarise(b1=median(b1), b2=median(b2), b3=median(b3), b4=median(b4)) %>%
    mutate(
        vis = factor(str_sub(visandsign, end=-9)),
        sign = factor(str_sub(visandsign, start=-8))
    )
pdfb = extract_samples(mcmcChain, ~ pred_y[i]) %>%
    join(pred_df, by="i")

#get typical mu values
typical_mu = extract_samples(mcmcChain, ~ typical_mu[visandsign], df) %>%
    spread(visandsign, typical_mu)

#differences between groups
typical_mu$group1 = rowMeans(typical_mu[,c("scatterplotnegative","scatterplotpositive","parallelCoordinatesnegative")])
typical_mu$group2 = rowMeans(typical_mu[,c("ordered_linepositive","donutnegative","stackedbarnegative","ordered_linenegative","stackedlinenegative","stackedareanegative")])
typical_mu$group3 = rowMeans(typical_mu[,c("parallelCoordinatespositive","radarpositive","linepositive")])
typical_mu$group4 = rowMeans(typical_mu[,c("donutpositive", "linenegative", "radarnegative", "stackedareapositive", "stackedbarpositive", "stackedlinepositive")])
tm_group_comp = rbind(
    data.frame(comparison="2-1", difference=with(typical_mu, group2 - group1)),
    data.frame(comparison="3-2", difference=with(typical_mu, group3 - group2)),
    data.frame(comparison="4-3", difference=with(typical_mu, group4 - group3))
)
openGraph(5,5)
ggplot(tm_group_comp,
        aes(x=comparison, y=difference/log(2))) + 
    geom_violin(linetype=0, fill="skyblue") + 
    geom_hline(yintercept=0, lty="dashed") +
    stat_summary(fun.data="median_hilow", alpha=.99) +
    coord_flip() 
saveGraph("output/typical_mu-high_precision_group", "pdf")


#tau
taudf = extract_samples(mcmcChain, ~ tau[visandsign], df)
ggplot(taudf,
        aes(x=visandsign, y=sqrt(1/tau))) + 
    geom_violin(linetype=0, fill="skyblue") + 
    stat_summary(fun.data="median_hilow") +
    coord_flip() 

utaudf = extract_samples(mcmcChain, ~ u_tau[visandsign], df)
openGraph(7,5)
ggplot(filter(utaudf, visandsign %in% c("scatterplotpositive","scatterplotnegative","parallelCoordinatesnegative")),
        aes(x=visandsign, y=sqrt(1/u_tau))) + 
    geom_violin(linetype=0, fill="skyblue") + 
    stat_summary(fun.data="median_hilow", alpha=.99) +
    coord_flip() 
saveGraph("output/u_tau-high_precision_group", "pdf")

#differences between scatter plot and parallel coordinates
utaudf %>%
    mutate(sd=sqrt(1/u_tau)) %>% 
    select(-u_tau) %>%
    spread(visandsign, sd) %>%
    mutate(ratio=scatterplotpositive / parallelCoordinatesnegative) %>%
    ggplot(aes(x=1, y=ratio)) + 
        geom_violin(linetype=0, fill="skyblue") + 
        stat_summary(fun.data="median_hilow", alpha=.99) +
        coord_flip()



#typical mu
tmdf = extract_samples(mcmcChain, ~ typical_mu[visandsign], df)
tmdf$visandsign_bymean = with(tmdf, reorder(visandsign, -typical_mu, mean))
ggplot(tmdf,
        aes(x=visandsign_bymean, y=typical_mu)) + 
        geom_violin(linetype=0, fill="skyblue") + 
        geom_hline(yintercept=log(.45), lty="dashed") +
        stat_summary(fun.data="median_hilow") +
#        ylim(-3.5, 0) +
        coord_flip() +
        theme_bw()

#generate pairwise comparisons of adjacent conditions
tmdf_comparison = filter(tmdf, visandsign != "radarnegative") %>%
    arrange(visandsign_bymean)
names(tmdf_comparison) = paste0(names(tmdf_comparison), "1")
tmdf_2 = filter(tmdf, visandsign != "scatterplotpositive") %>%
    arrange(visandsign_bymean)
names(tmdf_2) = paste0(names(tmdf_2), "2")
tmdf_comparison = cbind(tmdf_comparison, tmdf_2) %>%
    mutate(
        comparison = factor(paste0(visandsign_bymean1, " - ", visandsign_bymean2)),
        typical_mu_difference = typical_mu2 - typical_mu1
    )
tmdf_comparison$comparison = reorder(tmdf_comparison$comparison, as.numeric(tmdf_comparison$visandsign_bymean1))

ggplot(tmdf_comparison,
        aes(x=comparison, y=typical_mu_difference)) + 
    geom_violin(linetype=0, fill="skyblue") + 
    geom_hline(yintercept=0, lty="dashed") +
    stat_summary(fun.data="median_hilow") +
#    geom_segment(data=tmdf_intervals, mapping=aes(x=visandsign_bymean, xend=visandsign_bymean, y=tm_lower, yend=tm_upper), size=1.25) +
#    geom_point(data=tmdf_intervals, mapping=aes(x=visandsign_bymean, y=tm_mean), size=3, shape=3) +
#        ylim(-3.5, 0) +
    coord_flip()


    

#log-space model posteriors
ggplot(
        df,
#        df,
        aes(x=r, 
            y=log(jnd),
            color=not_censored,
            group=NA
#            color=p_chance_cutoff 
#            color=jnd > .4 | jnd > .95 - r
#            color=approach
            )) + 
    geom_point(alpha=.1) + 
    geom_hline(yintercept=log(.45), lty="dashed") + 
    stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
    stat_function(fun=function(x) log(x), lty="dashed", color="black") + 
    geom_abline(data=bdf_m, mapping=aes(intercept=b1, slope=b2)) +
#    stat_function(fun=function(x) log(.95 - x), lty="dashed", color="black") + 
    #    stat_smooth(method=rq, se=FALSE) +  
    stat_smooth(method=lm) + 
    facet_wrap(~visandsign)


#log-space prediction posteriors
ggplot(
        pdfb,
#        df,
        aes(x=r, 
            y=log(pred_y)
#            color=p_chance_cutoff 
#            color=jnd > .4 | jnd > .95 - r
#            color=approach
        )) + 
	geom_bin2d(breaks=list(
			x=seq(0.25, .85, by=.1), 
			y=seq(min(log(pdfb$pred_y)), max(log(pdfb$pred_y)), length=50)
		), mapping=aes(alpha=..density.., fill=1)) +
#    geom_point(alpha=.1) + 
    geom_hline(yintercept=log(.45), lty="dashed") + 
    geom_hline(yintercept=log(.40), lty="dashed") + 
    stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
    geom_abline(data=bdf_m, mapping=aes(intercept=b1, slope=b2), color="blue") +
#    stat_function(fun=function(x) log(.95 - x), lty="dashed", color="black") + 
    #    stat_smooth(method=rq, se=FALSE) +  
#    stat_smooth(method=lm) + 
    geom_point(data=df, aes(x=r, y=log(jnd)), alpha=.1) + 
    facet_wrap(~visandsign)




#log-space model posteriors
pred_r_step = .1/8
pred_r_min = .3
pred_r_max = .8

rdf = ldply(seq(pred_r_min, pred_r_max, pred_r_step), function(r) {
        within(bdf, {
            r <- r
            logy <- b1 + b2*r
        })
    }) %>%
    select(-b3, -b4) %>%	#don't need b3 / b4 for this
    mutate(
        vis = factor(str_sub(visandsign, end=-9)),
        sign = factor(str_sub(visandsign, start=-8))
        )
rdf_m = rdf %>% 
    group_by(visandsign, r) %>%
    summarise(logy=mean(logy))


ggplot(
        rdf,
#        df,
        aes(x=r, 
            y=logy,
#            color=p_chance_cutoff 
#            color=jnd > .4 | jnd > .95 - r
#            color=approach
        )) + 
	geom_bin2d(breaks=list(
			x=seq(pred_r_min - pred_r_step/2, pred_r_max + pred_r_step/2, pred_r_step), 
			y=seq(min(rdf$logy), max(rdf$logy), length=100)
		), mapping=aes(alpha=..density.., fill=1)) +
#    geom_point(alpha=.1) + 
    geom_hline(yintercept=log(.45), lty="dashed") + 
    geom_hline(yintercept=log(.40), lty="dashed") + 
    stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
    geom_abline(data=bdf_m, mapping=aes(intercept=b1, slope=b2), color="blue") +
#    stat_function(fun=function(x) log(.95 - x), lty="dashed", color="black") + 
    #    stat_smooth(method=rq, se=FALSE) +  
#    stat_smooth(method=lm) + 
    geom_point(data=df, mapping=aes(x=r, y=log(jnd), color=not_censored), alpha=.1) + 
    facet_wrap(~visandsign)


#log-space fit lines, single plot, no data
openGraph(9, 5.5)
ggplot(
        rdf,
        aes(x=r, 
            y=logy/log(2),
            group=visandsign
        )) + 
	geom_bin2d(breaks=list(
			x=seq(pred_r_min - pred_r_step/2, pred_r_max + pred_r_step/2, pred_r_step), 
			y=seq(min(rdf$logy/log(2)), max(rdf$logy/log(2)), length=150)
		), mapping=aes(alpha=..density.., fill=vis, color=NULL)) +
    geom_abline(data=bdf_m, mapping=aes(intercept=b1/log(2), slope=b2/log(2), color=vis, linetype=sign), size=1) + 
    geom_text(data=bdf_m, mapping=aes(y=I(b1/log(2) + .9 * b2/log(2)), x=.9, color=vis, label=visandsign)) +
    geom_hline(yintercept=log2(.45), lty="dashed") + 
    stat_function(fun=function(x) log2(1 - x), lty="dashed", color="black") +
    scale_alpha_continuous(range=c(0.01,0.6)) + 
    scale_y_continuous(labels=trans_format(function(x) 2^x, math_format(.x))) +
    scale_x_continuous(breaks=seq(0.3,0.8,by=0.1)) +
    annotation_logticks(sides="l")
saveGraph("output/final-model-log-space", "pdf")
    
    
#linear-space fit lines, single plot, no data
ggplot(
        rdf,
        aes(x=r, 
            y=exp(logy),
            color=visandsign
        )) + 
	geom_bin2d(breaks=list(
			x=seq(pred_r_min - pred_r_step/2, pred_r_max + pred_r_step/2, pred_r_step), 
			y=seq(exp(min(rdf$logy)), exp(max(rdf$logy)), length=200)
		), mapping=aes(alpha=..density.., fill=visandsign, color=NULL)) +
    scale_alpha_continuous(range=c(0.05,1)) + 
    geom_line(data=rdf_m, aes(x=r, y=exp(logy), color=visandsign), size=1) +  
    geom_hline(yintercept=.45, lty="dashed") + 
    geom_abline(slope=-1, intercept=1, lty="dashed") +
    ylim(0, .75)


#mean jnd (log space) for a particular r value
plot_r = .8
rdf$visandsign_bymean = with(rdf, reorder(visandsign, -logy, mean))

rdf_intervals = ddply(rdf[rdf$r==plot_r,], ~ visandsign_bymean, summarize, 
    logy_mean=mean(logy), logy_lower=quantile(logy,.025), logy_upper=quantile(logy,.975))

ggplot(rdf[rdf$r==plot_r,], aes(x=visandsign_bymean, y=logy)) + 
    geom_violin(linetype=0, fill="skyblue") + 
    geom_hline(yintercept=log(.45), lty="dashed") +
    geom_segment(data=rdf_intervals, mapping=aes(x=visandsign_bymean, xend=visandsign_bymean, y=logy_lower, yend=logy_upper), size=1.25) +
    geom_point(data=rdf_intervals, mapping=aes(x=visandsign_bymean, y=logy_mean), size=3, shape=3) +
    ylim(-3.5, 0) +
    coord_flip() +
    theme_bw()
grid.edit("geom_point.points", grep = TRUE, gp = gpar(lwd = 3))


#save
save.image(file=paste("output/censored_regression-random_effects-intercept-FINAL", 
				(if (final_model) "final" else "not_final"),
				".RData", sep=""))
#load("output/censored_regression-random_effects-intercept-FINALfinal.RData")

