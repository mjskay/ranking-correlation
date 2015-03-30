library(coda)
library(rjags)
library(runjags)
library(survival)
library(gamlss)
library(gamlss.cens)
library(ggplot2)
library(grid)
library(plyr)
library(dplyr)
library(stringr)
library(scales)

memory.limit(8000)

source("src/openGraphSaveGraph.R")

source("src/clean-data.R")

graphics.off()
fileNameRoot="output/censored_regression_nonrandom" # for constructing output filenames

#------------------------------------------------------------------------------
# THE MODEL.
include_predictions = TRUE
final_model = TRUE
include_typical = TRUE

modelstring = paste("
model {
	#MODEL
	#core model
	for (i in 1:n_data) {
		# latent variable log-linear model
		mu[i] <- b[visandsign_number[i],1] + b[visandsign_number[i],2]*r[i] + 
		         b[visandsign_number[i],3]*approach_value[i] + b[visandsign_number[i],4]*approach_value[i]*r[i] +
				 u[participant_number[i]]
    	y[i] ~ dlnorm(mu[i], tau[visandsign_number[i]])

		# right-censoring at this value of r
		y_is_censored[i] ~ dinterval(y[i], censoring_threshold[i])
  	}
    
	#participant random effects
	for (p in 1:n_participants) {
		u[p] ~ dnorm(0, u_tau[participant_visandsign_number[p]])
	}

	#priors
	for (v in 1:n_visandsigns) {
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
	for (v in 1:n_visandsigns) {
		typical_mu[v] <- b[v,1] + b[v,2]*typical_r
	}" else "", "

	#PREDICTIONS
	", if (include_predictions) "
	for (i in 1:n_pred) {
		pred_u[i] ~ dnorm(0, u_tau[pred_visandsign_number[i]])
		pred_mu[i] <- b[pred_visandsign_number[i],1] + b[pred_visandsign_number[i],2]*pred_r[i] +
					pred_u[i]
    	pred_y[i] ~ dlnorm(pred_mu[i], tau[pred_visandsign_number[i]])
	}" else "", "

}")
writeLines(modelstring,con="model.txt")


#------------------------------------------------------------------------------
# THE DATA.

#convert factors to 1-indexed numeric lists
participants = levels(factor(df$participant))
df$participant_number = as.numeric(factor(df$participant))
n_participants = max(df$participant_number)
visandsigns = levels(df$visandsign)
df$visandsign_number = as.numeric(df$visandsign)
n_visandsigns = max(df$visandsign_number)

#determine the visandsign each participant was assigned to
participant_visandsign_number = laply(1:n_participants, function(p) {
        filter(df, participant_number == p)[1,]$visandsign_number
    })

#right-censor the jnd values
y = df$jnd
y_is_censored = !df$not_censored
y[y_is_censored] = NA

#code approach as sum-to-zero
df$approach_value = ifelse(df$approach == "above", -1, 1)

#build data matrix
n_data = nrow(df)

dataList = list(
		n_data = n_data,
		y = y,
        y_is_censored = as.numeric(y_is_censored),
        r = df$r,
        approach_value = df$approach_value,
        visandsign_number = df$visandsign_number,
        n_visandsigns = n_visandsigns,
        censoring_threshold = df$censoring_threshold,
    	participant_number = df$participant_number,
    	n_participants = n_participants,
        participant_visandsign_number = participant_visandsign_number
    )


#------------------------------------------------------------------------------
# DATA FOR PREDICTIONS
if (include_predictions) {
    #predictions
    pred_df = expand.grid(
        visandsign_number=unique(df$visandsign_number), 
        r=unique(df$r))

    n_pred = nrow(pred_df)
	dataList$pred_visandsign_number = pred_df[,1]
	dataList$pred_r = pred_df[,2]
	dataList$n_pred = n_pred
}

#------------------------------------------------------------------------------
# INTIALIZE THE CHAINS.
#get maximum likelihood estimates of model parameters
coefs = ddply(df, ~ visandsign, function(df) {
        m = with(df, gamlss(Surv(censored_jnd, not_censored) ~ r * approach_value, family=cens(LOGNO)))
        c(
            sigma = exp(coef(m, "sigma")), 
            coef(m)
        )
    })

#y initialization should provide values for the latent variable only for observations 
#above the censoring threshold (values below the threshold are observed data so
#cannot be initialized here)
y_init = df$jnd
y_init[!y_is_censored] = NA

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
    jagsModel = run.jags("model.txt", data=dataList, monitor=parameters, initlist=inits_list, 
        method="parallel")
} else {
#    jagsModel = autorun.jags("model.txt", data=dataList, monitor=parameters, initlist=inits_list,
#        method="parallel", thin.sample=TRUE)    
    jagsModel = run.jags("model.txt", data=dataList, monitor=parameters, initlist=inits_list,
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
thinIndex = ceiling(seq(1, nrow(mcmcChain), by=1))#length=500))

# Extract chain values:
extract_2d_vector_sample = function(variable_name, n, n_name) {
	#extract a 2d vector sample in long-format where the row in the sample is the visandsign
	sample = NULL
	for (v in 1:n_visandsigns) {
		for (i in 1:n) {
			sdf = data.frame(v = mcmcChain[thinIndex, paste(variable_name, "[", v, ",", i,"]", sep="")])
			names(sdf) = variable_name
			sdf$visandsign = visandsigns[v]
			sdf[[n_name]] = i
			sample = rbind(sample, sdf)
		}
	}
	sample$visandsign = factor(sample$visandsign, levels=visandsigns)
	sample
}
extract_2d_vector_sample_wide = function(variable_name, n) {
	#extract a 2d vector sample in wide format where the row in the sample is the visandsign
	sample = NULL
	for (v in 1:n_visandsigns) {
        sdf = NULL
		for (i in 1:n) {
			sdf2 = data.frame(v = mcmcChain[thinIndex, paste(variable_name, "[", v, ",", i,"]", sep="")])
			names(sdf2) = paste0(variable_name, i)
            sdf = if (is.null(sdf)) sdf2 else cbind(sdf, sdf2)
		}
		sdf$visandsign = visandsigns[v]
		sample = rbind(sample, sdf)
	}
	sample$visandsign = factor(sample$visandsign, levels=visandsigns)
	sample
}
extract_1d_vector_sample = function(row_index, variable_name, n) {
	sample = NULL
	for (i in 1:n) {
		sdf = data.frame(v = mcmcChain[thinIndex, paste(variable_name, "[", i, "]", sep="")])
		names(sdf) = variable_name
        sdf = cbind(row_index[i,], sdf)
		sample = rbind(sample, sdf)
	}
	sample
}
extract_sample = function(mcmcChain, spec) {
    #first, extract the sample into a data frame
    variable_name = as.character(spec[[2]][[2]])
    long_df = extract_sample_long(mcmcChain,
        variable_name, 
        as.character(spec[[2]][-1:-2]))
    
    if (is.null(long_df$..)) {
        long_df
    } else {
        #if any of the columns were named "..", use it to form a wide version of the data
        wide_df = filter(long_df, .. == 1)
        wide_df$.. = NULL
        wide_df[[paste0(variable_name, "1")]] = wide_df[[variable_name]]
        wide_df[[variable_name]] = NULL
        for (i in 1:max(long_df$..)) {
            col_df = filter(long_df, .. == i)
            wide_df[[paste0(variable_name, i)]] = col_df[[variable_name]]
        }
        wide_df
    }   
}

extract_sample_long = function(mcmcChain, variable_name, index_names) {
    ldply(colnames(mcmcChain), function(colname) {
                    colname_parts = strsplit(colname,"\\,|\\]|\\[")[[1]]
                    if (colname_parts[1] == variable_name) {	#this is one of the variables we want 
                            indices = as.list(as.numeric(colname_parts[-1]))
                            names(indices) = index_names
                            values = list(mcmcChain[,colname])
                            names(values) = variable_name
                            data.frame(indices, values)
                    }            
            })
}
bdf = extract_sample(mcmcChain, ~ b[visandsign_number, ..]) %>%
    mutate(visandsign = factor(visandsigns[visandsign_number])) %>%
    select(-visandsign_number)
bdf_m = ddply(bdf, ~ visandsign, summarize, b1=median(b1), b2=median(b2), b3=median(b3)) %>%
    mutate(
        vis = factor(str_sub(visandsign, end=-9)),
        sign = factor(str_sub(visandsign, start=-8))
    )
pred_df$visandsign = factor(pred_df$visandsign_number, labels=visandsigns)
pdfb = extract_1d_vector_sample(pred_df, "pred_y", n_pred)

#get typical mu values
typical_mu = extract_sample(mcmcChain, ~ typical_mu[..])
for (i in 1:length(visandsigns)) {
    typical_mu[[visandsigns[[i]]]] = typical_mu[[paste0("typical_mu", i)]]
    typical_mu[[paste0("typical_mu", i)]] = NULL
}

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
saveGraph("output/mu-group-comparison", "pdf")


#tau
taudf=extract_sample(mcmcChain, ~ tau[visandsign_number])
taudf$visandsign = factor(visandsigns[taudf$visandsign_number])
ggplot(taudf,
        aes(x=visandsign, y=sqrt(1/tau))) + 
    geom_violin(linetype=0, fill="skyblue") + 
    stat_summary(fun.data="median_hilow") +
    coord_flip() 

utaudf=extract_sample(mcmcChain, ~ u_tau[visandsign_number])
utaudf$visandsign = factor(visandsigns[utaudf$visandsign_number])
openGraph(7,5)
ggplot(filter(utaudf, visandsign %in% c("scatterplotpositive","scatterplotnegative","parallelCoordinatesnegative")),
        aes(x=visandsign, y=sqrt(1/u_tau))) + 
    geom_violin(linetype=0, fill="skyblue") + 
    stat_summary(fun.data="median_hilow", alpha=.99) +
    coord_flip() 
saveGraph("output/u_tau-high_precision_group", "pdf")



#typical mu
tmdf=extract_sample(mcmcChain, ~ typical_mu[visandsign_number])
tmdf$visandsign = factor(visandsigns[tmdf$visandsign_number])
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
saveGraph("output/final-model-log-space.pdf")
    
    
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


b0_sample_full = NULL
for (j in 1:(n_acceptable_levels - 1)) {
	b0_sample_full = cbind(b0_sample_full, mcmcChain[, paste("b0[", j ,"]", sep="")])
}
b_survey_part_sample_full = NULL
s_survey_part_sample_full = NULL
for (j in 1:n_survey_parts) {
	b_survey_part_sample_full = cbind(b_survey_part_sample_full, mcmcChain[, paste("b_survey_part[", j ,"]", sep="")])
	s_survey_part_sample_full = cbind(s_survey_part_sample_full, mcmcChain[, paste("s_survey_part[", j ,"]", sep="")])
}
if (include_predictions) {
	pred_y_sample_full = NULL
	for (j in 1:n_pred) {
		pred_y_sample_full = cbind(pred_y_sample_full, mcmcChain[, paste("pred_y[", j ,"]", sep="")])
	}
}
b_sample_full = mcmcChain[, "b"] 
alpha_sample_full = mcmcChain[, "alpha"] 
model_number_sample = mcmcChain[, "model_number"]
p_model = laply(1:3, function(model_number) {
	sum(model_number_sample == model_number) / length(model_number_sample)
})

#subset to preferred model
best_model = which(p_model == max(p_model))
b0_sample = b0_sample_full[model_number_sample == best_model,]
b_survey_part_sample = b_survey_part_sample_full[model_number_sample == best_model,]
s_survey_part_sample = s_survey_part_sample_full[model_number_sample == best_model,]
b_sample = b_sample_full[model_number_sample == best_model]
alpha_sample = as.matrix(alpha_sample_full[model_number_sample == best_model])
if (include_predictions) pred_y_sample = pred_y_sample_full[model_number_sample == best_model,]
chain_length = nrow(b0_sample)

#save
save.image(file=paste("output/censored_regression-random_effects-intercept-FINAL", 
				(if (final_model) "final" else "not_final"),
				".RData", sep=""))
#load("output/censored_regression-random_effects-intercept-FINALfinal.RData")

