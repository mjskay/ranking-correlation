# Examine convergence of the Bayesian (censored, log-linear) model
# 
# Author: mjskay
###############################################################################

library(runjags)
library(coda)

load(file="output/bayesian_model.RData")

fit = add.summary(fit)

summary(fit)

plotfun = function(...) plot(fit, ..., layout=c(2,2), plot.type=c("trace", "ecdf", "histogram", "autocorr", "crosscorr"))
plotfun(vars="^b", file="output/model-params-b.pdf")
plotfun(vars="^tau", file="output/model-params-tau.pdf")
plotfun(vars="^u_tau", file="output/model-params-u_tau.pdf")
if (include_typical) plotfun(vars="^typ", file="output/model-params-typical_mu.pdf")
if (include_predictions) plotfun(vars="^pred", file="output/model-params-pred.pdf")
