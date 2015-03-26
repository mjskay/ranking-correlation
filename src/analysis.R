library(ggplot2)
library(plyr)
library(lme4)
library(visreg)
library(MASS)
library(quantreg)
library(gamlss)

#DATA CLEANING

#read data
df = read.csv("data/master.csv")

#some summary stats
visandsign_summaries = ddply(df, ~ visandsign, summarize, 
    p_chance = sum(jnd > .45)/length(jnd)	#proportion of jnd worse than chance
)
df = join(df, visandsign_summaries)
visandsign_rbase_summaries = ddply(df, ~ visandsign + rbase + approach, summarize, 
    jnd_mad = mad(jnd),
    jnd_median = median(jnd)
)
df = join(df, visandsign_rbase_summaries)

#adjusted r
ra_summaries = ddply(df, ~ visandsign + rbase, summarize, 
    jnd_mean_vs = mean(jnd),
    jnd_median_vs = median(jnd)
)
df = join(df, ra_summaries)
df$ra = df$rbase + ifelse(df$approach == "above", 0.5, -0.5) * df$jnd_median_vs  

#columns for filters used by original paper
df$p_chance_cutoff = df$p_chance > .2
df$mad_cutoff = abs(df$jnd - df$jnd_median) > 3 * df$jnd_mad 



#ANALYSIS 
#how much data was dropped due to the mad cutoff?
prop.table(xtabs(~ mad_cutoff + p_chance_cutoff, data=df), 2)



ggplot(
#    df[!df$mad_cutoff,], 
    df,
    aes(x=ordered(rbase), 
#        y=jnd, 
        y=log(jnd),
        fill=p_chance_cutoff)) + geom_violin() + facet_wrap(~visandsign)

# linear scale
ggplot(
        filter(df),
        aes(x=r, 
            y=jnd, 
#            color=p_chance_cutoff,
             color=jnd > .45,
#            color=jnd > .4 | jnd > .95 - rbase
            group=NA
    )) + 
    geom_point(size=1.5, alpha=.25) + 
#    geom_jitter(size=.75, position=position_jitter(width=.01, height=0)) + 
#    geom_line(alpha=.1, mapping=aes(group=participant:approach)) +
    geom_hline(yintercept=.45, lty="dashed") + 
    geom_abline(slope=-1, intercept=1, lty="dashed") +
    geom_abline(slope=1, intercept=0, lty="dashed") +
    xlim(0,1) +
    stat_smooth(method=lm) +
#    stat_smooth(method=gamlss, family=LNO) +   
#    stat_smooth(method=rq, se=FALSE) +  
#    stat_smooth(method=glm, family=Gamma(link=log)) +  
#    scale_y_log10() +
    facet_wrap(~visandsign)

#residual fit w.r.t. r (linear model) for one visualization



theme_set(theme_get() + theme(
    text=element_text(size=16, family="Avenir"),
    axis.text=element_text(size=rel(15/16)),
    axis.ticks.length=unit(8, "points"),
    line=element_line(size=1)
))


ggplot(
        filter(df, visandsign == "parallelCoordinatesnegative" & !mad_cutoff),
        aes(x=r, 
            y=jnd,
            color=approach
            )) + 
    geom_point(size=2, alpha=.25) +  
    stat_smooth(method=lm, se=FALSE) +
    stat_smooth(method=lm, se=FALSE, mapping=aes(group=NA))
     



#log scale
ggplot(
        df,
#        df,
        aes(x=r, 
            y=log(jnd),
            color=p_chance_cutoff
#            color=not_censored
#            color=jnd > .4 | jnd > .95 - rbase
#            color=approach
    )) + 
    geom_point(alpha=.1) + 
#    geom_line(alpha=.1, mapping=aes(group=participant:approach)) +
    geom_hline(yintercept=log(.45), lty="dashed") + 
    geom_hline(yintercept=log(.40), lty="dashed") + 
    stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
    stat_function(fun=function(x) log(x), lty="dashed", color="black") +
#    geom_ribbon(data=pdf2, mapping=aes(x=r, y=NULL, ymax=logy.fit + 1.96 * logy.se.fit, ymin=logy.fit - 1.96 * logy.se.fit), fill="blue", alpha=.25) + 
#    geom_line(data=pdf2, mapping=aes(x=r, y=logy.fit), color="green") + 
#    geom_line(data=pdf, mapping=aes(x=r, y=logy), color="red") + 
    #        geom_abline(data=bdf_m, mapping=aes(intercept=b1, slope=b2)) +
#    stat_function(fun=function(x) log(.95 - x), lty="dashed", color="black") + 
    #    stat_smooth(method=rq, se=FALSE) +  
    stat_smooth(method=lm) + 
    facet_wrap(~visandsign)
#df$visandsign = with(df, reorder(visandsign, jnd > .45, sum))


options(contrasts=c("contr.sum", "contr.poly"))
m = lm(jnd ~ ra * visandsign, data=df[!df$p_chance_cutoff,])
mg = glm(jnd ~ ra * visandsign, data=df[!df$p_chance_cutoff,], family=Gamma(link=log))
ml = lm(log(jnd) ~ ra * visandsign, data=df[!df$p_chance_cutoff,])

mr = rq(jnd ~ ra * visandsign, data=df[!df$p_chance_cutoff,])

m.spp = lmer(jnd ~ ra + (rbase|participant), data=df[!df$p_chance_cutoff & !df$mad_cutoff & df$visandsign == "scatterplotpositive",])
m.sppl = lmer(log(jnd) ~ ra + (rbase|participant), data=df[!df$p_chance_cutoff & df$visandsign == "scatterplotpositive",])
m.sppg = glmer(jnd ~ ra + (rbase|participant), data=df[!df$p_chance_cutoff & df$visandsign == "scatterplotpositive",], family=Gamma(link=log))

m.spp = lm(jnd ~ ra, data=df[!df$p_chance_cutoff & !df$mad_cutoff & df$visandsign == "scatterplotpositive",])
m.sppl = lm(log(jnd) ~ ra, data=df[!df$p_chance_cutoff & df$visandsign == "scatterplotpositive",])
m.sppg = glm(jnd ~ ra, data=df[!df$p_chance_cutoff & df$visandsign == "scatterplotpositive",], family=Gamma(link=log))


m.spp = lm(log(jnd) ~ ra, data=df[!df$p_chance_cutoff & !df$mad_cutoff & df$visandsign == "scatterplotpositive",])

m.spp = lm(log(jnd) ~ ra, data=df[!df$p_chance_cutoff & !df$mad_cutoff & df$visandsign == "scatterplotpositive",])


df_ag = ddply(df[!df$p_chance_cutoff & !df$mad_cutoff & df$visandsign == "scatterplotpositive",], ~ rbase + approach, function(df) {
        data.frame(
            ra = mean(df$ra),
            jnd = mean(df$jnd) 
            )
    })

m = lm(jnd ~ ra, data=df_ag)

