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

memory.limit(6096)

source("src/openGraphSaveGraph.R")

source("src/clean-data.R")

graphics.off()
fileNameRoot="output/censored_regression_nonrandom" # for constructing output filenames

#------------------------------------------------------------------------------
# THE MODEL.
include_predictions = TRUE
final_model = TRUE

modelstring = paste("
model {
	#core model
	for (i in 1:n_data) {
		# latent variable log-linear model
		mu[i] <- b[visandsign_number[i],1] + b[visandsign_number[i],2]*r[i] + 
		         b[visandsign_number[i],3]*approach_value[i] + b[visandsign_number[i],4]*approach_value[i]*r[i] 
    	y[i] ~ dlnorm(mu[i], tau[visandsign_number[i]])

		# right-censoring at this value of r
		y_is_censored[i] ~ dinterval(y[i], censoring_threshold[i])
  	}
    
	for (v in 1:n_visandsigns) {
		#priors for coefficents of latent model
		b[v,1] ~ dnorm(log(.45), 1E-2)	#prior on intercept is chance
		b[v,2] ~ dnorm(0, 1E-2)
		b[v,3] ~ dnorm(0, 1E-2)
		b[v,4] ~ dnorm(0, 1E-2)

    	#prior on variance of latent variable
    	tau[v] ~ dgamma(0.01, 0.01)
	}

	#typical performance depending on distribution of r
	typical_r ~ dunif(0.3,0.8)
	for (v in 1:n_visandsigns) {
		typical_mu[v] <- b[v,1] + b[v,2]*typical_r
	}

	", if (include_predictions) "
	for (i in 1:n_pred) {
		pred_mu[i] <- b[pred_visandsign_number[i],1] + b[pred_visandsign_number[i],2]*pred_r[i]
    	pred_y[i] ~ dlnorm(pred_mu[i], tau[pred_visandsign_number[i]])
	}" else "", "

}")
writeLines(modelstring,con="model.txt")


#------------------------------------------------------------------------------
# THE DATA.

#convert factors to 1-indexed numeric lists
visandsigns = levels(df$visandsign)
df$visandsign_number = as.numeric(df$visandsign)
n_visandsigns = max(df$visandsign_number)

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
        censoring_threshold = df$censoring_threshold
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
            tau = 1/exp(coef(m, "sigma"))^2, 
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
    tau=coefs$tau,
    y=y_init
)

#------------------------------------------------------------------------------
# RUN THE CHAINS

parameters = c("b" , "tau", "typical_r", "typical_mu")  # The parameter(s) to be monitored.
if (include_predictions) parameters = c(parameters, "pred_y") 


if (!final_model) {
    jagsModel = run.jags("model.txt", data=dataList, monitor=parameters, initlist=inits_list)
} else {
    jagsModel = autorun.jags("model.txt", data=dataList, monitor=parameters, initlist=inits_list,
        method="parallel", thin.sample=TRUE)    
}

codaSamples = as.mcmc.list(jagsModel)

# Create, initialize, and adapt the model:
#jagsModel = jags.model( "model.txt" , data=dataList, initlist=inits_list , 
#                        n.chains=nChains , n.adapt=adaptSteps )
#
## Burn-in:
#cat( "Burning in the MCMC chain...\n" )
#update( jagsModel , n.iter=burnInSteps )
#
## The saved MCMC chain:
#cat( "Sampling final MCMC chain...\n" )
#codaSamples = coda.samples( jagsModel , variable.names=parameters , 
#                            n.iter=nPerChain , thin=thinSteps )
# resulting codaSamples object has these indices: 
#   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]


#------------------------------------------------------------------------------
# EXAMINE THE RESULTS

checkConvergence = FALSE
if ( checkConvergence ) {
  show( summary( codaSamples ) )
  openGraph()
  plot( codaSamples , ask=F )  
  openGraph()
  autocorr.plot( codaSamples , ask=T )
}

# Convert coda-object codaSamples to matrix object for easier handling.
# But note that this concatenates the different chains into one long chain.
# Result is mcmcChain[ stepIdx , paramIdx ]
mcmcChain = as.matrix( codaSamples )
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
bdf_m = ddply(bdf, ~ visandsign, summarize, b1=median(b1), b2=median(b2), b3=median(b3))
pred_df$visandsign = factor(pred_df$visandsign_number, labels=visandsigns)
pdfb = extract_1d_vector_sample(pred_df, "pred_y", n_pred)

#get typical mu values
typical_mu = extract_sample(mcmcChain, ~ typical_mu[..])
for (i in 1:length(visandsigns)) {
    typical_mu[[visandsigns[[i]]]] = typical_mu[[paste0("typical_mu", i)]]
    typical_mu[[paste0("typical_mu", i)]] = NULL
}


tmdf=extract_sample(mcmcChain, ~ typical_mu[visandsign_number])
tmdf$visandsign = factor(visandsigns[tmdf$visandsign_number])

ggplot(extract_sample(mcmcChain, ~ typical_mu[i]),
    
    
    )


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
pred_r_min = .1
pred_r_max = .9

rdf = ldply(seq(pred_r_min, pred_r_max, pred_r_step), function(r) {
        within(bdf, {
            r <- r
            logy <- b1 + b2*r
        })
    })
rdf_m = ddply(rdf, ~ visandsign + r, summarize, logy=mean(logy))


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
ggplot(
        rdf,
        aes(x=r, 
            y=logy,
            color=visandsign
        )) + 
	geom_bin2d(breaks=list(
			x=seq(pred_r_min - pred_r_step/2, pred_r_max + pred_r_step/2, pred_r_step), 
			y=seq(min(rdf$logy), max(rdf$logy), length=150)
		), mapping=aes(alpha=..density.., fill=visandsign, color=NULL)) +
    scale_alpha_continuous(range=c(0.01,0.6)) + 
    geom_abline(data=bdf_m, mapping=aes(intercept=b1, slope=b2, color=visandsign), size=1) + 
    geom_hline(yintercept=log(.45), lty="dashed") + 
    stat_function(fun=function(x) log(1 - x), lty="dashed", color="black")
    
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

#get sample summaries
sampleSummary = summary(codaSamples)

#subset codaSamples to only those in the best model
codaSamplesBest = codaSamples
for (i in 1:length(codaSamplesBest)) {
	codaSamplesBest[[i]] = coda::as.mcmc(codaSamplesBest[[1]][codaSamplesBest[[1]][,"model_number"] == best_model,])
}
sampleSummaryBest = summary(codaSamplesBest)

#get predictions
if (include_predictions) {
	pred_df$pred_y_mean = laply(1:n_pred, function (i) {				
		sampleSummaryBest$statistics[paste("pred_y[", i, "]", sep=""),"Mean"]		
		})
	pred_df$pred_y_median = laply(1:n_pred, function (i) {				
		sampleSummaryBest$quantiles[paste("pred_y[", i, "]", sep=""),"50%"]
		})
	pred_df$pred_y_hdi = laply(1:n_pred, function (i) {				
		#highest-density interval
		hi = hist(pred_y_sample[,i], breaks=0:7 + 0.5)
		hi$mids[which(hi$density == max(hi$density))]
		})
	sqrt(mean((pred_df$acceptable_number - pred_df$pred_y_mean)^2, na.rm=TRUE))
	mean(abs(pred_df$acceptable_number - pred_df$pred_y_median), na.rm=TRUE)
	sqrt(mean((pred_df$acceptable_number - mean(pred_df$acceptable_number, na.rm=TRUE))^2, na.rm=TRUE))
	mean(abs(pred_df$acceptable_number - median(pred_df$acceptable_number, na.rm=TRUE)), na.rm=TRUE)
	
	write.csv(pred_df, file=paste("weather_rain_regression-scale_new-",
							(if (final_model) "final" else "not_final"), 
							"_predictions.csv", sep=""))
}


#save
save.image(file=paste("output/censored_regression-no_random_effects-", 
				(if (final_model) "final" else "not_final"),
				".RData", sep=""))
#load("output/censored_regression-no_random_effects-final.RData")

#source("weather_rain_regression_parameter_plots.R")
