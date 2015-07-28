# Misc plots of the Bayesian (censored, log-linear) model, most of which didn't make
# the paper (see the end of README.md for plots from the paper).
#
# No guarantees about readability here.
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
library(ggplot2)
library(grid)
library(stringr)
library(scales)


memory.limit(8000)


load(file="output/bayesian_model.RData")

#apply variable type prototypes from the origcinal data so that
#extract_samples will automatically convert variable indices based on factors
#back to named factors with the same levels as in the original data
fit %<>% apply_prototypes(df)

#extract samples of linear model coefficients
b_samples = extract_samples(fit, b[visandsign, ..]) %>%
    mutate( #split vis and sign apart so we can apply aesthetics to them separately
        vis = factor(str_sub(visandsign, end=-9)),
        sign = factor(str_sub(visandsign, start=-8))
    )
b_medians = b_samples %>%
    group_by(visandsign, vis, sign) %>%	#vis, sign is redundant here but we want to keep them around
    summarise(b1 = median(b1), b2 = median(b2), b3 = median(b3), b4 = median(b4))

#extract samples for variables that were indexed by visandsign
samples_by_visandsign = extract_samples(fit, cbind(tau, u_tau)[visandsign])

#differences in typical_mu between partialling-ranked groups of visandsigns
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
    ggeye(aes(x=comparison, y=difference/log(2))) + 
    geom_hline(y=0, lty="dashed")


#tau
samples_by_visandsign %>%
    ggeye(aes(x=visandsign, y=sqrt(1/tau)))

#participant tau for high-performing visandsigns
samples_by_visandsign %>%
    filter(visandsign %in% c("scatterplotpositive","scatterplotnegative","parallelCoordinatesnegative")) %>%
    ggeye(aes(x=visandsign, y=sqrt(1/u_tau)))

#differences between visandsigns
samples_by_visandsign %>% 
    filter(visandsign %in% c("scatterplotpositive","scatterplotnegative","parallelCoordinatesnegative")) %>%
    mutate(sd=sqrt(1/u_tau)) %>%
    compare_levels(sd, by=visandsign) %>% 
    ggeye(aes(x=visandsign, y=sd)) + 
    geom_hline(y=0, linetype="dashed")

#typical mu
typical_mu_samples = extract_samples(fit, typical_mu[visandsign]) %>%
    mutate(visandsign = reorder(visandsign, typical_mu, mean))
#estimated typical_mu in each visandsign
typical_mu_samples %>% 
    ggeye(aes(x=visandsign, y=typical_mu)) +  
    geom_hline(y=log(.45), lty="dashed")
#comparisons of typical_mu in each visandsign by order of increasing typical_mu
typical_mu_samples %>%
    compare_levels(typical_mu, by=visandsign, comparison=ordered) %>%
    ggeye(aes(x=visandsign, y=typical_mu)) +  
    geom_hline(y=0, lty="dashed")


#log-space fit lines, with data and censoring indicated
df %>%
    ggplot(aes(x=r, y=log(jnd),
            color=censored, group=NA)) + 
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
    geom_point(data=df, aes(x=r, y=log(jnd), color=censored), alpha=.1) + 
    facet_wrap(~visandsign)

#log-space fit lines with posterior density of mu, single plot, no data
mu_by_r %>%
    ggplot(aes(x=r, y=log_mu/log(2), group=visandsign)) + 
    geom_bin2d(breaks=list(
			x=seq(pred_r_min - pred_r_step/2, pred_r_max + pred_r_step/2, pred_r_step), 
			y=seq(min(mu_by_r$log_mu/log(2)), max(mu_by_r$log_mu/log(2)), length=150)
		), mapping=aes(alpha=..density.., fill=vis, color=NULL)) +
    geom_abline(data=b_medians, mapping=aes(intercept=b1/log(2), slope=b2/log(2), color=vis, linetype=sign), size=1) + 
    geom_text(data=b_medians, mapping=aes(y=I(b1/log(2) + .82 * b2/log(2)), x=.82, color=vis, label=visandsign), hjust=0) +
    geom_hline(yintercept=log2(.45), lty="dashed") + 
    scale_alpha_continuous(range=c(0.01,0.6), guide=FALSE) +
    scale_color_discrete(guide=FALSE) +
    scale_fill_discrete(guide=FALSE) +
    scale_y_continuous(labels=trans_format(function(x) 2^x, math_format(.x))) +
    scale_x_continuous(breaks=seq(0.3,0.8,by=0.1), limits=c(0.3,1)) +
    annotation_logticks(sides="l")

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
    ylim(0, .75)

#mean jnd (log space) at a particular r value for each visandsign
plot_r = .55
mu_by_r %>%
    mutate(visandsign = reorder(visandsign, log_mu, mean)) %>%
    filter(r == plot_r) %>%
    ggeye(aes(x=visandsign, y=log_mu)) + 
    geom_hline(y = log(.45), lty="dashed")

