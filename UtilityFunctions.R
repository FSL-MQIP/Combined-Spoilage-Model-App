##Utility functions
#Primary Growth Models 

Baranyi = function(LOG10N0, LOG10Nmax, lag, mumax, t) {
  ans <-  LOG10Nmax + log10((-1 + exp(mumax * lag) + exp(mumax * t))/(exp(mumax * t) - 1 + exp(mumax * lag) * 10^(LOG10Nmax - LOG10N0)))
  return(ans)
}

Buchanan = function(LOG10N0, LOG10Nmax, lag, mumax, t){
  ans <- LOG10N0 + (t >= lag) * (t <= (lag + (LOG10Nmax - LOG10N0) * log(10) / mumax)) * mumax * (t - lag) / log(10) + (t >= lag) * (t>(lag + (LOG10Nmax - LOG10N0) * log(10) / mumax)) * (LOG10Nmax - LOG10N0)
  return(ans)
}

Gompertz = function(LOG10N0, LOG10Nmax, lag, mumax, t){
  ans <- LOG10N0 + (LOG10Nmax - LOG10N0) * exp(-exp(mumax * exp(1) * (lag - t) / ((LOG10Nmax - LOG10N0) * log(10)) + 1))
  return(ans)
}


log10N <- function(LOG10N0, LOG10Nmax, lag, mumax, t, model_type) {
  if (model_type == "baranyi") {
    return(Baranyi(LOG10N0, LOG10Nmax, lag, mumax, t))
  }
  else if(model_type == 'buchanan') {
    return(Buchanan(LOG10N0, LOG10Nmax, lag, mumax, t))
  }
  else if(model_type == "gompertz"){
    return(Gompertz(LOG10N0, LOG10Nmax, lag, mumax, t))
  }else {
    stop(paste0(model_type, " is not a valid model type. Must be one of Buchanan, Baranyi, Gompertz"))
  }
}

#Secondary model 
muAtNewTemp <- function(newTemp, oldMu, oldTemp = 6, T0) { #T0 is -4.15 for ppc; T0 is -3.62 for spore
  numerator <- newTemp - T0
  denom <- oldTemp - T0
  newMu <- ((numerator / denom)^2) * oldMu
  
  return(newMu)
}

lagAtNewTemp <- function (newTemp, oldLag, oldTemp = 6, T0) { #T0 is -4.15 for ppc; T0 is -3.62 for spore
  numerator <- oldTemp -T0
  denom <- newTemp - T0
  newLag <- ( (numerator / denom)^2) * oldLag
  return(newLag)
}