# Back-of-the-napkin calculations to get some reasonable priors for the 
# Bayesian model based on Rensink & Baldridge.
#
# N.B. These calculations were used as a rough guide, the final priors
# and their rationale are discussed in the paper. 
# 
# Author: mjskay
###############################################################################

#derive some reasonable priors
Rensink_slope = -.24
Rensink_intercept = .24/.907
Rensink_jnd = function(r) Rensink_intercept + Rensink_slope * r 

#First, get approx slope and intercept in log space
Rensink_log_slope = (log(Rensink_jnd(.8)) - log(Rensink_jnd(.3))) / .5
Rensink_log_intercept =  log(Rensink_jnd(.3)) - Rensink_log_slope * .3
Rensink_log_jnd = function(r) Rensink_log_intercept + Rensink_log_slope * r

#let's get a prior on the variance. Assuming individual differences are large-ish
#(say, the sd is equal to abs(slope), such that it's reasonable to see a JND at r==.3 that
#is comparable to the mean JND at r==.8)
sd_1 = log(Rensink_jnd(.3)) - log(Rensink_jnd(.8))
#now let's assume the actual sd can be as high as 2 times this
sd_2 = sd_1 * 2
#prior on variance
tau_inverse_max = sd_2^2
tau_max = 1/tau_inverse_max

#Now a prior on the intercept: chance, with variance on the order of twice the distance to the scatteplot intercept
b_1 = log(.45)
b_1_sd = 2 * (log(.45) - Rensink_log_intercept)
b_1_tau = 1/(b_1_sd^2)

#prior on the slope: up to twice it
b_2 = 0
b_2_sd = 2 * abs(Rensink_log_slope)
b_2_tau = 1/(b_2_sd^2)

#prior on approach effect: approximately the same as Rensink at r = .3, up to 2 times it
b_3 = (log(.2) - log(.16))
b_3_sd = 2 * b_3
b_3_tau = 1/(b_3_sd^2)

#prior on approach interaction: 0, with same scale as approach effect
b_4 = 0
b_4_tau = b_3_tau


print(data.frame(
        tau_inverse_max,
        b_1, b_1_tau,
        b_2, b_2_tau,
        b_3, b_3_tau,
        b_4, b_4_tau
        ))
