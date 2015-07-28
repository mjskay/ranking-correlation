# Helper functions for plotting the fit of the (non-Bayesian) models fit by gamlss
# 
# Author: mjskay
###############################################################################

library(dplyr)
library(ggplot2)        # ggplot, geom_..., etc
library(scales)         # trans_format
library(grid)
library(gtable)


#model fit plotting function
plot_model_residuals = function(m, log_x=FALSE, log_base=2, data_color = "#999999", sd_color="#d62d0e", dnorm_color="#066fff") {
    
    data = data.frame(
        eval(m$call$data),
        `Fitted JND`=fitted(m), 
        `Normalized Quantile Residuals`=resid(m), 
        check.names=FALSE)
    
    scale_location = data %>%
        group_by(interval=cut(`Fitted JND`, 15)) %>%
        summarize(
            `Fitted JND` = mean(as.numeric(strsplit(as.character(interval), '(\\(|\\,|\\])')[[1]][2:3])),
            standard_deviation = sd(`Normalized Quantile Residuals`)
        )
    
    if (log_x) {
        #change log base for easier to understand better output
        data$`Fitted JND` = data$`Fitted JND` / log(log_base) 
        scale_location$`Fitted JND` = scale_location$`Fitted JND` / log(log_base)
    }
    
    fitted_plot = ggplot(data, 
            aes(x=`Fitted JND`, y=`Normalized Quantile Residuals`)
        ) +
        geom_point(color=data_color, size=2, alpha=0.25) +
        geom_rug(sides="r", color=data_color) +
        geom_hline(y=0, linetype="dashed", alpha=0.5) + 
        geom_linerange(data=scale_location, aes(ymin=standard_deviation, ymax=-standard_deviation, y=NULL), color=sd_color, size=1.0) +
        ylim(-5.5,5.5)
    if (log_x) {
        fitted_plot = fitted_plot +
            scale_x_continuous(
                limits=c(min(data$`Fitted JND`), max(data$`Fitted JND`)),
                labels=trans_format(function(x) log_base^x, math_format(.x))) +
            geom_vline(y=log2(.45)) +
            annotation_logticks(sides="b")
    } else {
        fitted_plot = fitted_plot +
            xlim(min(data$`Fitted JND`), max(data$`Fitted JND`))
    }
    fitted_plot = ggplotGrob(fitted_plot)
    
    
    residual_hist = ggplot(data, aes(x=`Normalized Quantile Residuals`)) +
        stat_density(fill=data_color) +
        geom_vline(x=0, linetype="dashed", alpha=0.5) + 
        xlim(-5.5,5.5) +
        ylim(0,0.6) +
        stat_function(fun=dnorm, linetype="dashed", color=dnorm_color) +
        coord_flip()
    residual_hist = ggplotGrob(residual_hist)
    
    grid.newpage()
    gtable(heights=unit(1,"null"), widths=unit(c(3,5,1),c("lines","null","null"))) %>%
        gtable_add_grob(fitted_plot[,1:3], 1, 1) %>%
        gtable_add_grob(fitted_plot[,4:5], 1, 2) %>%
        gtable_add_grob(residual_hist[,4:5], 1, 3) %>%
        grid.draw()
    
}

#model fitting plotting function for a given viz, showing residuals relative to r
plot_model_residuals_by_r = function(df, m, plotted_visandsign, log_y=FALSE, log_base=2, data_color = "#999999", sd_color = "#d62d0e") {
    #calculate residual variance
    residual_df = df %>%
        filter(visandsign == plotted_visandsign & !mad_cutoff) %>%
        cbind(prediction=predict(m, newdata=.)) %>%
        mutate(
            jnd = if (log_y) log(jnd, base=log_base) else jnd,
            prediction = if (log_y) prediction / log(log_base) else prediction
        ) %>%
        group_by(r) %>%
        summarize(
            sd = sd(jnd - prediction),
            jnd = mean(prediction)
        )
    
    #calculate fit line
    coefs = as.list(coef(m))
    fit_slope = coefs$r + coefs[[paste0("r:visandsign", plotted_visandsign)]]
    fit_intercept = coefs$`(Intercept)` + coefs[[paste0("visandsign", plotted_visandsign)]]
    
    #chance threshold
    chance = if (log_y) log(.45, base=log_base) else .45
    
    if (log_y) {
        df$jnd = log(df$jnd, base=log_base)
        fit_slope = fit_slope / log(log_base)
        fit_intercept = fit_intercept / log(log_base)
    }
    
    p = ggplot(
            filter(df, visandsign == plotted_visandsign & !mad_cutoff),
            aes(x=r, 
                y=jnd 
            )) + 
        geom_point(size=3, alpha=.25, color=data_color) + 
        geom_abline(slope=fit_slope, intercept=fit_intercept, size=1.0) +
        geom_linerange(data=residual_df, aes(ymin=jnd - sd, ymax=jnd + sd, y=NULL), color=sd_color, size=1.0)
    
    if (log_y) {
        p = p + 
            annotation_logticks(sides="l")
    }
    
    p
}
