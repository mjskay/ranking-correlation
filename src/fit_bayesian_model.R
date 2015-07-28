# Fit the Bayesian (censored, log-linear) model for JND ~ r
# 
# Author: mjskay
###############################################################################

library(plyr)           #**ply
library(dplyr)          #filter, group_by, mutate, select, etc
library(coda)           #autocorr.plot
library(runjags)        #run.jags
library(tidybayes)      #compose_data, apply_prototypes, extract_samples, compare_levels
                        #see https://github.com/mjskay/tidybayes
library(metabayes)      #metajags
                        #see https://github.com/mjskay/metabayes


memory.limit(8000)


# MODEL SPECIFICATION
# flags to allow smaller models to be built for testing purposes. Set to FALSE for faster runs
include_predictions = TRUE     #should we include predictions in the model? 
final_model = TRUE             #should we run the final (long) model?
include_typical = TRUE         #should we calculate typical mean JND (in log space) assuming some distribution of r?

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


# LOAD THE DATA
source("src/clean_data.R")

#determine the visandsign each participant was assigned to
#N.B. this relies on the assumption (correct for the Harrison data)
#that all participants were assigned to only one visandsign
participant_visandsign = df %>%
    group_by(participant) %>%
    slice(1) %>%
    with(as.numeric(visandsign))

#build data list for input into sampler
data_list = df %>%
    mutate(y = ifelse(censored, NA, jnd)) %>%	#values of latent variable are observed data only if not censored
    select(y, r, approach_value, visandsign, participant, censoring_threshold, censored) %>%
    compose_data(participant_visandsign)


# DATA FOR PREDICTIONS
if (include_predictions) {
    #add values of independent variables we want to predict from to the data list
    pred_df = expand.grid(
        pred_visandsign=1:nlevels(df$visandsign), 
        pred_r=unique(df$r))
    data_list = data_list %>% 
        compose_data(pred=pred_df)
    #restructure pred_df so that its columns are useful later (when bound onto predicted data)
    pred_df = pred_df %>%
        mutate(i=1:nrow(.), visandsign=factor(pred_visandsign, labels=levels(df$visandsign))) %>%    
        select(r=pred_r, visandsign, i)
}


# FIT THE MODEL
#determine the parameters to be monitored in the fit
parameters = c("b", "tau", "u_tau", 
    if (include_predictions) "pred_y", 
    if (include_typical) c("typical_r", "typical_mu")
)

fit = if (!final_model) {
        run.jags(code(model), data=data_list, monitor=parameters, 
            method="parallel")
    } else {
        run.jags(code(model), data=data_list, monitor=parameters,
    	    adapt=5000, burnin = 100000, sample = 10000, thin=10,
            summarise = FALSE,  #don't keep summary data in fit (makes saved model smaller and can be generated later)
            method="parallel")
    }


#SAVE MODEL
save.image(file=paste0(
        "output/bayesian_model", 
        (if (!include_predictions) "_~predictions"),
		(if (!final_model) "_~final"),
        (if (!include_typical) "_~typical"),        
		".RData"))
