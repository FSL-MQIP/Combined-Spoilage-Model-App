# Load packages
library(dplyr)
library(tidyr)
library(fitdistrplus)
library(tibble)
library(EnvStats)         # to load rtri function 
library(truncnorm)        # to load rtruncnorm function
library(jmuOutlier)       # to load rlaplace function
library(formula.tools)    # to load 'rhs'
library(purrr)            # to load 'map'
library(deSolve)          # to load 'ode'
library(rlang)

# Load utility functions
source("UtilityFunctions_dynamic_growth.R")

# Input data 
# Initial contamination levels
ppc_initial <- read.csv("InputData/ppc_initial.csv")
spore_initial <- read.csv("InputData/spore_initial.csv")
no_spoil_initial <- read.csv("InputData/initialcounts_060622.csv")
# Growth parameters
data_ppc <- read.csv("InputData/ppc_gp.csv")
data_sporeformer <- read.csv("InputData/sporeformer_gp.csv")
gp_noSpoil <- read.csv("InputData/NoSpoil_gp.csv")
colnames(gp_noSpoil)[1] <- "isolate"
# AT frequency
ppc_ST_freq <- read.csv("InputData/ppc_STfreq_clean.csv")
spore_AT_freq <- read.csv("InputData/spore_ATfreq_clean.csv")

# Set seed
set.seed(1)

# Generate input parameter distributions 
# Generate uncertainty distributions for initial contamination levels
# PPC initial contamination by bootstrapping 
means <- numeric(1000)
std_deviations <- numeric(1000)

for (i in 1:1000) {
  sample <- sample(ppc_initial$LOG.initial, size = nrow(ppc_initial), replace = TRUE)
  sample_mean <- mean(sample)
  sample_std <- sd(sample)
  means[i] <- sample_mean
  std_deviations[i] <- sample_std
}

ppc_mean_nfit  <- fitdist(means, "norm")
ppc_sd_nfit  <- fitdist(std_deviations, "norm")
mean_ppc_mean <- ppc_mean_nfit$estimate[1]
sd_ppc_mean <- ppc_mean_nfit$estimate[2]
mean_ppc_sd <- ppc_sd_nfit$estimate[1]
sd_ppc_sd <- ppc_sd_nfit$estimate[2]

# Generate normal distribution for ppc_mean 
ppc_mean_distr <- rnorm(1000, mean_ppc_mean, sd_ppc_mean)
# Generate normal distribution for ppc_sd
ppc_sd_distr <- rnorm(1000, mean_ppc_sd, sd_ppc_sd)

# Sporeformer initial contamination by bootstrapping
means <- numeric(1000)
std_deviations <- numeric(1000)

for (i in 1:1000) {
  sample <- sample(spore_initial$log10MPN, size = nrow(spore_initial), replace = TRUE)
  sample_mean <- mean(sample)
  sample_std <- sd(sample)
  means[i] <- sample_mean
  std_deviations[i] <- sample_std
}

spore_mean_nfit  <- fitdist(means, "norm")
spore_sd_nfit  <- fitdist(std_deviations, "norm")
mean_spore_mean <- spore_mean_nfit$estimate[1]
sd_spore_mean <- spore_mean_nfit$estimate[2]
mean_spore_sd <- spore_sd_nfit$estimate[1]
sd_spore_sd <- spore_sd_nfit$estimate[2]

# Generate normal distribution for spore_mean 
spore_mean_distr <- rnorm(1000, mean_spore_mean, sd_spore_mean)
# Generate normal distribution for spore_sd
spore_sd_distr <- rnorm(1000, mean_spore_sd, sd_spore_sd)

# Get mean and sd for PPC growth parameters h0, b, log10Nmax and Topt
data_ppc$h0 <- data_ppc$lag * data_ppc$mumax
Tmin_PPC <- -4.15
T_PPC <- 6
data_ppc$b <- sqrt(data_ppc$mumax/log(10)*24)/(T_PPC - Tmin_PPC) # mumax in log10/day
result_df <- data_ppc %>%
  group_by(isolate) %>%
  summarise(
    Mean_h0 = mean(h0),
    StdDev_h0 = sd(h0),
    Mean_b = mean(b),
    StdDev_b = sd(b),
    Mean_LOG10Nmax = mean(LOG10Nmax),
    StdDev_LOG10Nmax = sd(LOG10Nmax),
    Topt = Topt
  ) %>%
  distinct()
result_df <- as.data.frame(result_df)
result_df$STorAT <- unique(data_ppc$STorAT)
result_df$STorAT <- paste0("ST_", result_df$STorAT)

# Get mean and sd for sporeformer growth parameters h0, b, log10Nmax and Topt
data_sporeformer$h0 <- data_sporeformer$lag * data_sporeformer$mumax
Tmin_spore <- -3.62
T_spore <- 6
data_sporeformer$b <- sqrt(data_sporeformer$mumax/log(10))/(T_spore - Tmin_spore) # mumax in log10/day
result_df_1 <- data_sporeformer %>%
  group_by(isolate) %>%
  summarise(
    Mean_h0 = h0,
    StdDev_h0 = mean(result_df$StdDev_h0),
    Mean_b = b,
    StdDev_b = mean(result_df$StdDev_b),
    Mean_LOG10Nmax = LOG10Nmax,
    StdDev_LOG10Nmax = mean(result_df$StdDev_LOG10Nmax),
    Topt = Topt
  ) %>%
  distinct()
result_df_1 <- as.data.frame(result_df_1)
result_df_1$STorAT <- unique(data_sporeformer$STorAT)
result_df_1$STorAT <- paste0("AT_", result_df_1$STorAT)

# Join growth parameter data 
growth_parameter <- rbind(result_df,result_df_1,gp_noSpoil)

# Set up dataframe for simulation (100 lots, 100 units)
n_sim <- 100
n_unit <- 100
lot_id <- rep(seq(1, n_sim), each = n_unit)
unit_id <- rep(seq(1,n_unit), times = n_sim)
data <- data.frame(lot_id, unit_id)

# Temperature profile
# Stage 1: facility storage 
## (a)  Sample the temperature distribution
data$T_F <- rep(runif(n_sim*n_unit,min=3.5,max=4.5)) #uniform distribution
## (b) Sample the storage time (in days) distribution
data$t_F <- rep(runif(n_sim*n_unit,min=1,max=2)) #uniform distribution

# Stage 2: transport from facility to retail store
## (a)  Sample the temperature distribution
data$T_T <- rep(rtri(n_sim*n_unit,min=1.7,max=10.0,mode=4.4)) #triangular distribution
## (b) Sample the transportation time (in days) distribution
data$t_T <- rep(rtri(n_sim*n_unit,min=1,max=10,mode=5))

# Stage 3: storage/display at retail store
## (a)  Sample the temperature distribution
data$T_S <- rep(rtruncnorm(n_sim*n_unit,a=-1.4,b=5.4,mean=2.3,sd=1.8)) #truncated normal distribution
## (b) Sample the storage time (in days) distribution
data$t_S <- rep(rtruncnorm(n_sim*n_unit,a=0.042,b=10.0, mean=1.821,sd=3.3)) #truncated normal distribution

## Stage 4: transportation from retail store to home
## (a)  Sample the temperature distribution
data$T_T2 <- rep(rtruncnorm(n_sim*n_unit,a=0,b=10,mean=8.5,sd=1.0)) #truncated normal distribution
## (b) Sample the transportation time (in days) distribution 
data$t_T2 <- rep(rtruncnorm(n_sim*n_unit,a=0.01,b=0.24, mean=0.04,sd=0.02)) #truncated normal distribution

## Stage 5: home storage 
## (a)  Sample the temperature distribution
temps <- rep(NA, 10000)
for (i in 1:10000){
  number <- rlaplace(1,m=4.06,s=2.31)
  while (number > 15 | number < -1) {
    number <- rlaplace(1,m=4.06,s=2.31) #truncated laplace distribution 
  }
  temps[i] <- number
}
data$T_H <- temps
## (b) Define t_H as 35 days for all units
data$t_H <- rep(35, each = n_sim*n_unit)

# Generate spoilage frequency and assign spoilage types
num_ppc <- round(40/ 100 * n_sim * n_unit)
num_spore <- round(40/ 100 * n_sim * n_unit)
num_no_spoil <- n_sim*n_unit - num_ppc - num_spore
spoiler_types <- c(rep("PPC", num_ppc), rep("Spore", num_spore), rep("No Spoil", num_no_spoil))
data$spoilage_type <- sample(spoiler_types, n_sim * n_unit)

# Assign AT/ST types 
# ppc
model_data_ppc <- subset(data, spoilage_type == "PPC")
model_data_ppc$STorAT <- NA
for (i in 1:num_ppc) {
  sampled_value <- sample(ppc_ST_freq$ClosestST, 1, replace = T)
  while (sampled_value == "9_23") {
    sampled_value <- sample(ppc_ST_freq$ClosestST, 1, replace = T)
  }
  model_data_ppc$STorAT[i] <- sampled_value
}
model_data_ppc$STorAT <- paste0("ST_", model_data_ppc$STorAT)
# spore
model_data_spore <- subset(data, spoilage_type == "Spore")
model_data_spore$STorAT <- sample(spore_AT_freq$ClosestAT, num_spore,replace=T)
model_data_spore$STorAT <- paste0("AT_", model_data_spore$STorAT)
# no spoil
model_data_noSpoil <- subset(data, spoilage_type == "No Spoil")
model_data_noSpoil$STorAT <- rep("NS_0", num_no_spoil)

# Assign initial contamination distributions
# ppc
ppc_initial_mean <- sample(ppc_mean_distr, 100)
ppc_initial_sd <- sample(ppc_sd_distr, 100)
df1 <- data.frame(
  lot_id = unique(model_data_ppc$lot_id),
  initial_mean = rep(ppc_initial_mean, length.out = length(unique(model_data_ppc$lot_id))),
  initial_sd = rep(ppc_initial_sd, length.out = length(unique(model_data_ppc$lot_id)))
)
model_data_ppc <- merge(model_data_ppc, df1, by = "lot_id")
model_data_ppc <- model_data_ppc %>%
  rowwise() %>%
  mutate(N0 = rnorm(n = 1, mean = initial_mean, sd = initial_sd))
model_data_ppc <- as.data.frame(model_data_ppc)

# sporeformer 
spore_initial_mean <- sample(spore_mean_distr, 100)
spore_initial_sd <- sample(spore_sd_distr, 100)
df2 <- data.frame(
  lot_id = unique(model_data_spore$lot_id),
  initial_mean = rep(spore_initial_mean, length.out = length(unique(model_data_spore$lot_id))),
  initial_sd = rep(spore_initial_sd, length.out = length(unique(model_data_spore$lot_id)))
)
model_data_spore <- merge(model_data_spore, df2, by = "lot_id")
model_data_spore <- model_data_spore %>%
  rowwise() %>%
  mutate(N0 = rnorm(n = 1, mean = initial_mean, sd = initial_sd))
model_data_spore <- as.data.frame(model_data_spore)

# no Spoil
no_spoil_initial <- no_spoil_initial %>%
  filter(Day_Initial_Actual<=7)%>%
  filter(spoilagetype_actual=="noSpoil")
noSpoil_nfit  <- fitdist(log10(no_spoil_initial$SPC_DI), "norm")
model_data_noSpoil$initial_mean <- noSpoil_nfit$estimate[1]
model_data_noSpoil$initial_sd <- noSpoil_nfit$estimate[2]
model_data_noSpoil <- model_data_noSpoil %>%
  rowwise() %>%
  mutate(N0 = rnorm(n = 1, mean = initial_mean, sd = initial_sd))
model_data_noSpoil <- as.data.frame(model_data_noSpoil)

# join data
model_data <- rbind(model_data_ppc,model_data_spore,model_data_noSpoil)

# Generate allele index
model_data$allele_index <- match(model_data$STorAT, growth_parameter$STorAT)

# Assign growth parameters 
model_data$Mean_h0 <- growth_parameter$Mean_h0[model_data$allele_index]
model_data$StdDev_h0 <- growth_parameter$StdDev_h0[model_data$allele_index]
model_data$Mean_b <- growth_parameter$Mean_b[model_data$allele_index]
model_data$StdDev_b <- growth_parameter$StdDev_b[model_data$allele_index]
model_data$Mean_Nmax <- growth_parameter$Mean_LOG10Nmax[model_data$allele_index]
model_data$StdDev_Nmax <- growth_parameter$StdDev_LOG10Nmax[model_data$allele_index]
model_data$Topt <- growth_parameter$Topt[model_data$allele_index]

model_data <- model_data %>%
  rowwise() %>%
  mutate(h0 = rtruncnorm(n = 1,
                         a = 0,
                         mean = Mean_h0, 
                         sd = StdDev_h0))

model_data$Q0 <- 1/(exp(model_data$h0)-1)

model_data <- model_data %>%
  rowwise() %>%
  mutate(b = rtruncnorm(n = 1,
                        a = 0,
                        mean = Mean_b, 
                        sd = StdDev_b))

model_data <- model_data %>%
  rowwise() %>%
  mutate(Nmax = rnorm(n = 1,
                         mean = Mean_Nmax, 
                         sd = StdDev_Nmax))

model_data <- model_data %>%
  mutate(Tmin = case_when(
    spoilage_type == "Spore" ~ -3.62,
    spoilage_type == "PPC" ~ -4.15,
    spoilage_type == "No Spoil" ~ 0
  ))

model_data$mu_opt <- (model_data$b*(model_data$Topt-model_data$Tmin))^2 

model_data <- model_data %>%
  mutate(mu_opt = ifelse(spoilage_type == "No Spoil", 0, mu_opt),
         Q0 = ifelse(spoilage_type  == "No Spoil", 0, Q0),
         b = ifelse(spoilage_type == "No Spoil", 0, b))

# Model temperature profiles of 10000 units HTST milk 
env_cond_time <- matrix(c(rep(0,10000),
                          model_data$t_F, 
                          model_data$t_F+0.001,
                          model_data$t_F + model_data$t_T,
                          model_data$t_F + model_data$t_T+0.001,
                          model_data$t_F + model_data$t_T + model_data$t_S,
                          model_data$t_F + model_data$t_T + model_data$t_S+0.001,
                          model_data$t_F + model_data$t_T + model_data$t_S + model_data$t_T2,
                          model_data$t_F + model_data$t_T + model_data$t_S + model_data$t_T2+0.001,
                          model_data$t_F + model_data$t_T + model_data$t_S + model_data$t_T2 + model_data$t_H), ncol = 10)

env_cond_temp <- matrix(c(model_data$T_F, 
                          model_data$T_F,
                          model_data$T_T,
                          model_data$T_T,
                          model_data$T_S,
                          model_data$T_S,
                          model_data$T_T2,
                          model_data$T_T2,
                          model_data$T_H,
                          model_data$T_H), ncol = 10)

# Run simulation
final_conc <- model_data %>%
  rowwise() %>%
  mutate(final_conc_isolate = {
    my_primary <- list(mu_opt = mu_opt, Nmax = Nmax, N0 = N0, Q0 = Q0)
    sec_temperature <- list(model = "reducedRatkowsky", xmin = Tmin, b = b, xopt = Topt)
    my_secondary <- list(temperature = sec_temperature)
    growth <- predict_dynamic_growth(times = env_cond_time,
                                     env_conditions = tibble(time = env_cond_time,
                                                             temperature = env_cond_temp),
                                     my_primary,
                                     my_secondary)
    sim <- growth$simulation
    return(tail(sim$logN, 1))
  }) %>%
  pull(final_conc_isolate)
