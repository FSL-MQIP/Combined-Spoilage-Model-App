## Utility Functions
# Last edited: 03/15/23

## (1) baranyi_log10N
## (2) baranyi_without_lag_log10N
# Source: Equations for the 2 models below are copied from the Nlms package in R (https://rdrr.io/cran/nlsMicrobio/src/R/growthmodels.R)

# Purpose: Calculate log10N using Baranyi model 

# Parameter for baranyi_log10N
# (i) t: time in hours
# (ii) lag: length of lag phase
# (iii) mumax: growth rate (ln/hour)
# (iv) LOG10N0: initial microbial concentration 
# (v) LOG10Nmax: carrying capacity

# Parameters for baranyi_without_lag
# (i) t: time in hours
# (ii) mumax: growth rate (ln/hour)
# (iii) LOG10N0: initial microbial concentration 
# (iv) LOG10Nmax: carrying capacity

# Functions: 
baranyi_log10N = function(t,lag,mumax,LOG10N0,LOG10Nmax) {
  ans <- LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0)))
  return(ans)
}


baranyi_without_lag_log10N = function(t,mumax,LOG10N0,LOG10Nmax) {
  ans <- (LOG10Nmax - log10(1 + (10^(LOG10Nmax - LOG10N0) - 1) * exp(-mumax * t)))
  return(ans)
}
