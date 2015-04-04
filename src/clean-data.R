### DATA CLEANING

library(dplyr)


#read data
df = read.csv("data/master.csv")

df = df %>%
    #clean up columns from original data
    rename(
        r = rbase
    ) %>%

    #calculate values used to derive filtering conditions used in original paper
    group_by(visandsign) %>%
    mutate(
        p_chance = sum(jnd > .45)/length(jnd),	#proportion of jnd worse than chance
        p_chance_cutoff = p_chance > .2			#visandsigns with > 20% observations of jnd worse than chance
    ) %>%
    group_by(visandsign, r, approach) %>%
    mutate(
        #observations > 3 median-absolute deviations from the median within each group
        mad_cutoff = abs(jnd - median(jnd)) > 3 * mad(jnd)
    ) %>%
    ungroup() %>%

    #censorship based on being above / below ceiling/floor, and chance
    mutate(
        censoring_threshold = ifelse(approach == "below", 
            pmin(r - .05, .4), 
            pmin(.95 - r, .4)),
        not_censored = jnd <= censoring_threshold,
        censored_jnd = pmin(jnd, censoring_threshold)
    )

#approach should be coded as sum-to-zero so that other coefficients can be 
#interpreted as relative to the mean of both approaches. We relevel first 
#so that "below" is assigned 1 and "above" -1 (this is just to eliminate
#a double-negative so that the sign of the coefficient of approach is positive,
#makes interpretation slightly simpler)
df$approach = relevel(df$approach, "below")
contrasts(df$approach) = contr.sum
