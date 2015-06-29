library(coda)           #autocorr.plot
library(runjags)        #run.jags
library(ggplot2)        #ggplot, stat_..., geom_..., etc
library(plyr)           #**ply
library(dplyr)          #filter, group_by, mutate, select, etc
library(tidybayes)      #compose_data, apply_prototypes, extract_samples, compare_levels
library(metabayes)      #metajags
library(stringr)        #str_sub

memory.limit(8000)

source("src/openGraphSaveGraph.R")

source("src/clean-data.R")

#------------------------------------------------------------------------------
# THE MODEL.
include_predictions = FALSE
final_model = FALSE
include_typical = FALSE

model = metajags({
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
    if (include_typical) {
    	typical_r ~ dunif(0.3,0.8)
    	for (v in 1:n_visandsign) {
    		typical_mu[v] <- b[v,1] + b[v,2]*typical_r
    	}
    }

	#PREDICTIONS
    if (include_predictions) {
    	for (i in 1:n_pred) {
    		pred_u[i] ~ dnorm(0, u_tau[pred_visandsign[i]])
    		pred_mu[i] <- b[pred_visandsign[i],1] + b[pred_visandsign[i],2]*pred_r[i] +
    					pred_u[i]
        	pred_y[i] ~ dlnorm(pred_mu[i], tau[pred_visandsign[i]])
    	}
    }
})

#------------------------------------------------------------------------------
# THE DATA.

#determine the visandsign each participant was assigned to
#N.B. this relies on the assumption (correct for the Harrison data)
#that all participants were assigned to only one visandsign
participant_visandsign = df %>%
    group_by(participant) %>%
    slice(1) %>%
    with(as.numeric(visandsign))

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
        pred_visandsign=1:nlevels(df$visandsign), 
        pred_r=unique(df$r))
    data_list = data_list %>% 
        compose_data(pred=pred_df)
    #make pred_df more useful for later
    pred_df = pred_df %>%
        mutate(i=1:nrow(.), visandsign=factor(pred_visandsign, labels=levels(df$visandsign))) %>%    
        select(r=pred_r, visandsign, i)
}

#------------------------------------------------------------------------------
# RUN THE CHAINS

#determine the parameters to be monitored in the fit
parameters = c("b", "tau", "u_tau", 
    if (include_predictions) "pred_y", 
    if (include_typical) c("typical_r", "typical_mu")
)

#fit the model
fit = if (!final_model) {
        run.jags(code(model), data=data_list, monitor=parameters, 
            method="parallel")
    } else {
    #    autorun.jags(code(model), data=data_list, monitor=parameters,
    #        method="parallel", thin.sample=TRUE)    
        run.jags(code(model), data=data_list, monitor=parameters,
    	    adapt=5000, burnin = 100000, sample = 10000, thin=10,
            method="parallel")
    }


#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

#apply variable type prototypes from the original data so that
#extract_samples will automatically convert variable indices based on factors
#back to named factors with the same levels as in the original data
fit = fit %>%
    apply_prototypes(df)

checkConvergence = FALSE
if ( checkConvergence ) {
    plot(fit, vars="^b", file="output/model-params-b.pdf")
    plot(fit, vars="^tau", file="output/model-params-tau.pdf")
    plot(fit, vars="^u_tau", file="output/model-params-u_tau.pdf")
    plot(fit, vars="^typ", file="output/model-params-typical_mu.pdf")
    plot(fit, vars="^pred", file="output/model-params-pred.pdf")

    summary(fit)
    pdf(file="output/model-autocorr.pdf")
    autocorr.plot(as.mcmc.list(codaSamples), ask=FALSE)
    dev.off()
}

#extract samples of linear model coefficients
b_samples = extract_samples(fit, b[visandsign, ..]) %>%
    mutate( #split vis and sign apart so we can apply aesthetics to them separately
        vis = factor(str_sub(visandsign, end=-9)),
        sign = factor(str_sub(visandsign, start=-8))
    )
b_medians = b_samples %>%
    group_by(visandsign, vis, sign) %>%	#vis, sign is redundant here but we want to keep them around
    summarise(b1 = median(b1), b2 = median(b2), b3 = median(b3), b4 = median(b4)) %>%

#extract samples for variables that were indexed by visandsign
samples_by_visandsign = extract_samples(fit, cbind(tau, u_tau)[visandsign])

#differences in typical_mu between partialling-ranked groups of visandsigns
openGraph(5,5)
extract_samples(fit, typical_mu[visandsign] | visandsign) %>%
    mutate(     #group by partial ranking
        group1 = rowMeans(cbind(scatterplotnegative, scatterplotpositive, parallelCoordinatesnegative)),
        group2 = rowMeans(cbind(ordered_linepositive, donutnegative, stackedbarnegative, ordered_linenegative, stackedlinenegative, stackedareanegative)),
        group3 = rowMeans(cbind(parallelCoordinatespositive, radarpositive, linepositive)),
        group4 = rowMeans(cbind(donutpositive, linenegative, radarnegative, stackedareapositive, stackedbarpositive, stackedlinepositive))
    ) %>%
    with(rbind( #get differences between groups
        data.frame(comparison = "2-1", difference = group2 - group1),
        data.frame(comparison = "3-2", difference = group3 - group2),
        data.frame(comparison = "4-3", difference = group4 - group3)
    )) %>%
    ggraindrop(aes(x=comparison, y=difference/log(2))) + 
        geom_hline(y=0, lty="dashed")
saveGraph("output/typical_mu-group_differences", "pdf")


#tau
samples_by_visandsign %>%
    ggraindrop(aes(x=visandsign, y=sqrt(1/tau)))

#participant tau for high-performing visandsigns
openGraph(7,5)
samples_by_visandsign %>%
    filter(visandsign %in% c("scatterplotpositive","scatterplotnegative","parallelCoordinatesnegative")) %>%
    ggraindrop(aes(x=visandsign, y=sqrt(1/u_tau)))
saveGraph("output/u_tau-high_precision_group", "pdf")
#differences between visandsigns
samples_by_visandsign %>% 
    filter(visandsign %in% c("scatterplotpositive","scatterplotnegative","parallelCoordinatesnegative")) %>%
    mutate(sd=sqrt(1/u_tau)) %>%
    compare_levels(sd, by=visandsign) %>% 
    ggraindrop(aes(x=visandsign, y=sd)) + 
        geom_hline(y=0, linetype="dashed")

#typical mu
typical_mu_samples = extract_samples(fit, typical_mu[visandsign]) %>%
    mutate(visandsign = reorder(visandsign, typical_mu, mean))
#estimated typical_mu in each visandsign
typical_mu_samples %>% 
    ggraindrop(aes(x=visandsign, y=typical_mu)) +  
        geom_hline(y=log(.45), lty="dashed")
#comparisons of typical_mu in each visandsign by order of increasing typical_mu
typical_mu_samples %>%
    compare_levels(typical_mu, by=visandsign, comparison=ordered) %>%
    ggraindrop(aes(x=visandsign, y=typical_mu)) +  
        geom_hline(y=0, lty="dashed")

    
#log-space fit lines, with data and censoring indicated
df %>%
    ggplot(aes(x=r, y=log(jnd),
            color=not_censored, group=NA)) + 
        geom_point(alpha=.1) + 
        geom_hline(yintercept=log(.45), lty="dashed") + 
        stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
        stat_function(fun=function(x) log(x), lty="dashed", color="black") + 
        geom_abline(data=b_medians, mapping=aes(intercept=b1, slope=b2)) +
        facet_wrap(~visandsign)

#log-space posterior predictions
predictions = extract_samples(fit, pred_y[i]) %>%
    left_join(pred_df, by="i")
predictions %>%
    ggplot(aes(x=r, y=log(pred_y))) + 
    	geom_bin2d(breaks=list(
    			x=seq(0.25, .85, by=.1), 
    			y=seq(min(log(predictions$pred_y)), max(log(predictions$pred_y)), length=50)
    		), mapping=aes(alpha=..density.., fill=1)) +
        geom_point(data=df, aes(x=r, y=log(jnd)), alpha=.1) + 
        geom_hline(yintercept=log(.45), lty="dashed") +  
        stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
        stat_function(fun=function(x) log(x), lty="dashed", color="black") + 
        geom_abline(data=b_medians, mapping=aes(intercept=b1, slope=b2), color="red") +
        facet_wrap(~visandsign)

# generate predicted log(mu) for a set of values of r
pred_r_min = .3
pred_r_max = .8
pred_r_step = .1/8
mu_by_r = ldply(seq(pred_r_min, pred_r_max, pred_r_step), function(r) {
        within(b_samples, {
            r <- r
            log_mu <- b1 + b2*r
        })
    }) %>%
    select(-b3, -b4)    #don't need b3 / b4 for this
mu_median_by_r = mu_by_r %>% 
    group_by(visandsign, r) %>%
    summarise(log_mu=mean(log_mu))

#log-space fit lines with posterior density of mu, facetted, with data and censoring indicated
mu_by_r %>%
    ggplot(aes(x=r, y=log_mu)) + 
    	geom_bin2d(breaks=list(
			x=seq(pred_r_min - pred_r_step/2, pred_r_max + pred_r_step/2, pred_r_step), 
			y=seq(min(mu_by_r$log_mu), max(mu_by_r$log_mu), length=100)
		), mapping=aes(alpha=..density.., fill=1)) +
        geom_hline(yintercept=log(.45), lty="dashed") + 
        stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
        stat_function(fun=function(x) log(x), lty="dashed", color="black") + 
        geom_abline(data=b_medians, mapping=aes(intercept=b1, slope=b2), color="blue") +
        geom_point(data=df, aes(x=r, y=log(jnd), color=not_censored), alpha=.1) + 
        facet_wrap(~visandsign)

#log-space fit lines with posterior density of mu, single plot, no data
openGraph(9, 5.5)
mu_by_r %>%
    ggplot(aes(x=r, y=log_mu/log(2), group=visandsign)) + 
    	geom_bin2d(breaks=list(
			x=seq(pred_r_min - pred_r_step/2, pred_r_max + pred_r_step/2, pred_r_step), 
			y=seq(min(mu_by_r$log_mu/log(2)), max(mu_by_r$log_mu/log(2)), length=150)
		), mapping=aes(alpha=..density.., fill=vis, color=NULL)) +
        geom_abline(data=b_medians, mapping=aes(intercept=b1/log(2), slope=b2/log(2), color=vis, linetype=sign), size=1) + 
        geom_text(data=b_medians, mapping=aes(y=I(b1/log(2) + .9 * b2/log(2)), x=.9, color=vis, label=visandsign)) +
        geom_hline(yintercept=log2(.45), lty="dashed") + 
        scale_alpha_continuous(range=c(0.01,0.6)) + 
        scale_y_continuous(labels=trans_format(function(x) 2^x, math_format(.x))) +
        scale_x_continuous(breaks=seq(0.3,0.8,by=0.1)) +
        annotation_logticks(sides="l")
saveGraph("output/final-model-log-space", "pdf")
        
#linear-space fit lines with posterior density of mu, single plot, no data
mu_by_r %>%
    ggplot(aes(x=r, y=exp(log_mu), color=visandsign)) + 
    	geom_bin2d(breaks=list(
			x=seq(pred_r_min - pred_r_step/2, pred_r_max + pred_r_step/2, pred_r_step), 
			y=seq(exp(min(mu_by_r$log_mu)), exp(max(mu_by_r$log_mu)), length=200)
		), mapping=aes(alpha=..density.., fill=visandsign, color=NULL)) +
        scale_alpha_continuous(range=c(0.05,1)) + 
        geom_line(data=mu_median_by_r, aes(x=r, y=exp(log_mu), color=visandsign), size=1) +  
        geom_hline(yintercept=.45, lty="dashed") + 
        geom_abline(slope=-1, intercept=1, lty="dashed") +
        ylim(0, .75)

#mean jnd (log space) at a particular r value for each visandsign
plot_r = .55
mu_by_r %>%
    mutate(visandsign = reorder(visandsign, -log_mu, mean)) %>%
    filter(r == plot_r) %>%
    ggraindrop(aes(x=visandsign, y=log_mu)) + 
        geom_hline(y = log(.45), lty="dashed")

    
    
#save
save.image(file=paste("output/censored_regression-random_effects-intercept-FINAL", 
				(if (final_model) "final" else "not_final"),
				".RData", sep=""))
#load("output/censored_regression-random_effects-intercept-FINALfinal.RData")

