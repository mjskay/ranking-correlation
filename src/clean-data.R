library(ggplot2)
library(dplyr)
library(lme4)
library(visreg)
library(MASS)
library(gamlss)

#DATA CLEANING

#read data
df = read.csv("data/master.csv")

#clean up columns from original data
df = df %>%
    rename(r = rbase)

#calculate values used to derive filtering conditions used in original paper
df = df %>%
    group_by(visandsign) %>%
    summarize(p_chance = sum(jnd > .45)/length(jnd)) %>%	#proportion of jnd worse than chance
    join(df)

df = df %>%
    group_by(visandsign, r, approach) %>%
    summarize(
        jnd_mad = mad(jnd),			#within-group median absolute deviation 
        jnd_median = median(jnd)
    ) %>%
    join(df)

#filters used by original paper: when either is true, that data was excluded
df$p_chance_cutoff = df$p_chance > .2	#visandsigns with > 20% observations of jnd worse than chance
df$mad_cutoff = abs(df$jnd - df$jnd_median) > 3 * df$jnd_mad	#observations > 3 mads from the median

#adjusted r (TODO: ensure properly replicating old method)
df = df %>%
    group_by(visandsign, r) %>%
    summarize(jnd_mean_within_vsr = mean(jnd)) %>%
    join(df) %>%
    mutate(
        ra = r + ifelse(approach == "above", 0.5, -0.5) * jnd_mean_within_vsr
    )

#censorship based on being above / below ceiling/floor , and chance
df = mutate(df,
    censoring_threshold = ifelse(approach == "below", 
                    pmin(r - .05, .4), 
                    pmin(.95 - r, .4)),
    not_censored = jnd < censoring_threshold,
    censored_jnd = pmin(jnd, censoring_threshold)
)

#use sum-to-zero contrasts for approach so that we can estimate jnd as the 
#average of the above/below values
contrasts(df$approach) = contr.sum








