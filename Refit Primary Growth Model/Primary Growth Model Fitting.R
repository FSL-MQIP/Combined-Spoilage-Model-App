# Load packages 
library(nlsMicrobio)
library(minpack.lm)

# Import data 
# Sporeformers
data1 = read.csv("SporeformerRawGrowth.csv")
data1 = na.omit(data1)
colnames(data1)[1] = "t"
data2 = read.csv("SporeformerRawGrowthTrial2.csv")
data2 = na.omit(data2)
colnames(data2)[1] = "t"

# PPC
data3 = read.csv("PPCRawGrowthData.csv")
data3 = na.omit(data3)
colnames(data3)[1] = "t"
data3$LOG10N <- as.numeric(data3$LOG10N)

# Load primary growth models
source("UtilityFunctions_baranyi.R")

# Subset data
R103286Rep3 <- subset(data3, Isolate == "R10-3286" & Rep =="3")
plot(R103286Rep3$t,R103286Rep3$LOG10N)
mod <- lm(R103286Rep3$LOG10N[4:7] ~ R103286Rep3$t[4:7])
slope <- coef(mod) [2]
slope*2.303

# Fit Baranyi 
R103286Rep3.bar_LM <- nlsLM(LOG10N ~ baranyi_log10N(t,lag,mumax,LOG10N0,LOG10Nmax), data=R103286Rep3,
                      start=list (
                      LOG10N0 = 4,
                      lag = 12,
                      mumax = 0.07931, 
                      LOG10Nmax = 8),
                      lower = c(0,0,0,0))
                           
# Fitting Baranyi without lag 
H80237.bar_nolag_LM<- nlsLM(LOG10N ~ baranyi_without_lag_log10N(t,mumax,LOG10N0,LOG10Nmax), data=H80237,
                            start=list (
                            LOG10N0 = 3.5,
                            mumax = 0.3704, 
                            LOG10Nmax = 7.5), 
                            lower = c(0,0,0))




