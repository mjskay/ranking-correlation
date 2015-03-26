


#censored regression models
m.linear = gamlss(jnd ~ r * visandsign * approach,
    sigma.formula = ~ visandsign,
    data=df)

df2 = df[df$visandsign %in% c("ordered_linenegative","radarpositive" ),]
df2$visandsign = factor(df2$visandsign)

m.linear = gamlss(jnd ~ r * visandsign, sigma.formula = ~ visandsign, data=df2)

m.ln = gamlss(jnd ~ r * visandsign + approach, data=df[,c("jnd","r","visandsign","approach")], family=LOGNO)
m.ln.c = gamlss(Surv(censored_jnd, not_censored) ~ I(r - 0.3) * visandsign + approach_value,
    sigma.formula = ~ visandsign, 
    data=dfv, 
    family=cens(LOGNO))
m.ln.c.2 = survreg(Surv(censored_jnd, not_censored) ~ r * visandsign * approach_value, 
    data=dfv, 
    dist="lognormal")
m.n.c = gamlss(Surv(censored_jnd, not_censored) ~ r * visandsign + approach,
    sigma.formula = ~ visandsign, 
    data=df[c("censored_jnd","r","visandsign","approach","not_censored","participant")], 
    family=cens(NO))

pdf = expand.grid(r=seq(.3,.8,by=.1/8), approach_value=0, visandsign=levels(df$visandsign)) %>%
    cbind(logy=predict(m.ln.c, newdata=.))

pdf2 = expand.grid(r=seq(.3,.8,by=.1/8), approach_value=0, visandsign=levels(df$visandsign)) %>%
    cbind(logy=predict(m.ln.c.2, newdata=., se.fit=TRUE, type="link"))


#
mls = m.ln.c.2 %>%
    ref.grid(list(
            visandsign=levels(df$visandsign), 
            r=.8, 
            approach_value=0
            )) %>% 
    lsmeans(pairwise ~ visandsign) 


#plot, log scale, with data, faceted
ggplot(
        df,
#        df,
        aes(x=rbase, 
            y=log(jnd)
#            color=not_censored
#            color=jnd > .4 | jnd > .95 - rbase
#            color=approach
        )) + 
    geom_point(alpha=.1, mapping=aes(color=not_censored)) + 
    geom_hline(yintercept=log(.45), lty="dashed") + 
    geom_hline(yintercept=log(.40), lty="dashed") + 
    stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
    stat_function(fun=function(x) log(x), lty="dashed", color="black") +
    geom_ribbon(data=pdf2, mapping=aes(x=r, y=NULL, ymax=logy.fit + 1.96 * logy.se.fit, ymin=logy.fit - 1.96 * logy.se.fit), fill="blue", alpha=.25) + 
    geom_line(data=pdf2, mapping=aes(x=r, y=logy.fit), color="green") + 
#    geom_line(data=pdf, mapping=aes(x=r, y=logy), color="red") + 
    #        geom_abline(data=bdf_m, mapping=aes(intercept=b1, slope=b2)) +
#    stat_function(fun=function(x) log(.95 - x), lty="dashed", color="black") + 
    #    stat_smooth(method=rq, se=FALSE) +  
#    stat_smooth(method=lm) + 
    facet_wrap(~visandsign)
#df$visandsign = with(df, reorder(visandsign, jnd > .45, sum))


#plot, log scale, no data or facets
ggplot(
        pdf2,
#        df,
        aes(x=r, 
            y=logy.fit
#            color=not_censored
#            color=jnd > .4 | jnd > .95 - rbase
#            color=approach
        )) + 
    geom_hline(yintercept=log(.45), lty="dashed") + 
    geom_hline(yintercept=log(.40), lty="dashed") + 
    stat_function(fun=function(x) log(1 - x), lty="dashed", color="black") + 
    stat_function(fun=function(x) log(x), lty="dashed", color="black") +
    geom_ribbon(data=pdf2, mapping=aes(x=r, y=NULL, ymax=logy.fit + 1.96 * logy.se.fit, ymin=logy.fit - 1.96 * logy.se.fit, fill=visandsign), alpha=.25) + 
    geom_line(data=pdf2, mapping=aes(x=r, y=logy.fit, color=visandsign))
#    geom_line(data=pdf, mapping=aes(x=r, y=logy), color="red") + 
#        geom_abline(data=bdf_m, mapping=aes(intercept=b1, slope=b2)) +
#    stat_function(fun=function(x) log(.95 - x), lty="dashed", color="black") + 
#    stat_smooth(method=rq, se=FALSE) +  
#    stat_smooth(method=lm) + 
#df$visandsign = with(df, reorder(visandsign, jnd > .45, sum))




r. = .3
approach. = "below"
c_threshold = function(r) { 
    ifelse(approach. == "below", 
        pmin(r, .45), 
        pmin(1 - r, .45))
}
c_threshold_2 = function(r) { 
    ifelse(approach. == "below", 
        pmin(r - 0.05, .4), 
        pmin(.95 - r, .4))
}
ggplot(filter(df, r == r. & approach == approach.),
    aes(x=log(jnd), fill=p_chance_cutoff)
    ) +
    stat_bin() +
    geom_vline(x=log(c_threshold(r.))) +
#    geom_vline(x=log(c_threshold(r.)) - .2) +
    geom_vline(x=log(c_threshold_2(r.))) +
    facet_wrap(~visandsign)
    

plot(density(filter(df, r==.3, visandsign=="radarnegative")$jnd))