library(dplyr)

# Process VSL data 
combined_data <- read.csv("combined_data.csv")
combined_data$SampleID <- gsub("-\\d+", "", combined_data$SampleID)
# Filter out flavored milk samples 
filtered_combined_data <- subset(combined_data, !(endsWith(SampleID, "5") | endsWith(SampleID, "6")))
write.csv(filtered_combined_data,"filtered_combined_data.csv")

combined_data_sub <- subset(filtered_combined_data, Day_Initial <= 3)

combined_data_PPC <- subset(combined_data_sub,spoilagetype == "PPC")
CVTA_D7 <- log10(combined_data_PPC$CVTA_D7)
CVTA_D7 <- na.omit(CVTA_D7)
summary_stats_CVTA_D7 <- quantile(CVTA_D7, probs = c(0.25,0.5,0.75))
percentage_above_6_CVTA_D7 <- mean(CVTA_D7 > 6)*100

CVTA_D14 <- log10(combined_data_PPC$CVTA_D14)
CVTA_D14 <- na.omit(CVTA_D14)
summary_stats_CVTA_D14 <- quantile(CVTA_D14, probs = c(0.25,0.5,0.75))
percentage_above_6_CVTA_D14 <- mean(CVTA_D14 > 6)*100

CVTA_D21 <- log10(combined_data_PPC$CVTA_D21)
CVTA_D21 <- na.omit(CVTA_D21)
summary_stats_CVTA_D21 <- quantile(CVTA_D21, probs = c(0.25,0.5,0.75))
percentage_above_6_CVTA_D21 <- mean(CVTA_D21 > 6)*100

combined_data_Spore <- subset(combined_data_sub, spoilagetype %in% c("spore spoilage", "no spoilage"))
SPC_D1 <- combined_data_sub$SPC_DI
SPC_D7 <- log10(combined_data_Spore$SPC_D7)
SPC_D7 <- na.omit(SPC_D7)
summary_stats_SPC_D7 <- quantile(SPC_D7, probs = c(0.25,0.5,0.75))
percentage_above_6_SPC_D7 <- mean(SPC_D7 > 6)*100

SPC_D14 <- log10(combined_data_Spore$SPC_D14)
SPC_D14 <- na.omit(SPC_D14)
summary_stats_SPC_D14 <- quantile(SPC_D14, probs = c(0.25,0.5,0.75))
percentage_above_6_SPC_D14 <- mean(SPC_D14 > 6)*100

SPC_D21 <- log10(combined_data_Spore$SPC_D21)
SPC_D21 <- na.omit(SPC_D21)
summary_stats_SPC_D21 <- quantile(SPC_D21, probs = c(0.25,0.5,0.75))
percentage_above_6_SPC_D21 <- mean(SPC_D21 > 6)*100

# Validation 
# Load packages
library(tidyverse)
library(ggplot2)
library(gridExtra)
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
# Growth parameters
data_ppc <- read.csv("InputData/ppc_gp.csv")
data_ppc$mumax <- 0.684 * data_ppc$mumax
data_sporeformer <- read.csv("InputData/sporeformer_gp.csv")
# AT frequency
ppc_ST_freq <- read.csv("InputData/ppc_STfreq_clean.csv")
spore_AT_freq <- read.csv("InputData/spore_ATfreq_clean.csv")

# Set the seed for reproducibility
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
growth_parameter <- rbind(result_df,result_df_1)

# Set up dataframe for simulation (100 lots, 100 units)
n_sim <- 100
n_unit <- 100
lot_id <- rep(seq(1, n_sim), each = n_unit)
unit_id <- rep(seq(1,n_unit), times = n_sim)
data <- data.frame(lot_id, unit_id)

# Sanity check 
# data$P_ppc = 1  # validation for PPC
data$P_ppc = 0  # validation for Spore

# Assign spoilage type
result_list <- list()
for (i in 1:n_sim) {
  data_lot <- subset(data, lot_id == as.character(i))
  num_ppc <- round(unique(data_lot$P_ppc) * n_unit)
  num_spore_unspoil <- n_unit - num_ppc
  spoiler_types <- c(rep("PPC", num_ppc), rep("Spore_Unspoil", num_spore_unspoil))
  data_lot$spoilage_type <- sample(spoiler_types, n_unit)
  result_list[[i]] <- data_lot
}
data <- do.call(rbind,result_list)

# Assign AT/ST types 
# ppc
model_data_ppc <- subset(data, spoilage_type == "PPC")
model_data_ppc$STorAT <- NA
for (i in 1:nrow(model_data_ppc)) {
  sampled_value <- sample(ppc_ST_freq$ClosestST, 1, replace = T)
  while (sampled_value == "9_23") {
    sampled_value <- sample(ppc_ST_freq$ClosestST, 1, replace = T)
  }
  model_data_ppc$STorAT[i] <- sampled_value
}
model_data_ppc$STorAT <- paste0("ST_", model_data_ppc$STorAT)
# # spore
model_data_spore <- subset(data, spoilage_type == "Spore_Unspoil")
model_data_spore$STorAT <- sample(spore_AT_freq$ClosestAT, nrow(model_data_spore),replace=T)
model_data_spore$STorAT <- paste0("AT_", model_data_spore$STorAT)

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

# filter partial cells
model_data_ppc$N0_MPN_ppc <- 10^model_data_ppc$N0*1900
model_data_ppc_filtered <- subset(model_data_ppc, N0_MPN_ppc < 1)
model_data_ppc_filtered$new_column <- runif(nrow(model_data_ppc_filtered))
model_data_ppc_filtered$N0_MPN_ppc_assigned <- ifelse(model_data_ppc_filtered$N0_MPN_ppc > model_data_ppc_filtered$new_column, 1, 
                                                      ifelse(model_data_ppc_filtered$N0_MPN_ppc < model_data_ppc_filtered$new_column, 0, NA))
model_data_ppc_filtered <- model_data_ppc_filtered[, !(names(model_data_ppc_filtered) %in% c("N0_MPN_ppc", "new_column"))]
names(model_data_ppc_filtered)[names(model_data_ppc_filtered) == "N0_MPN_ppc_assigned"] <- "N0_MPN_ppc"
filtered_rows <- anti_join(model_data_ppc, model_data_ppc_filtered, by = c("lot_id", "unit_id"))
model_data_ppc <- rbind(filtered_rows,model_data_ppc_filtered)
model_data_ppc$N0 = model_data_ppc$N0_MPN_ppc/1900
model_data_ppc <- subset(model_data_ppc, select = -N0_MPN_ppc)
model_data_ppc <- as.data.frame(model_data_ppc)

# spore 
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

# filter partial spores
model_data_spore$N0_MPN_spore <- 10^model_data_spore$N0*1900
model_data_spore_filtered <- subset(model_data_spore, N0_MPN_spore < 1)
model_data_spore_filtered$new_column <- runif(nrow(model_data_spore_filtered))
model_data_spore_filtered$N0_MPN_spore_assigned <- ifelse(model_data_spore_filtered$N0_MPN_spore > model_data_spore_filtered$new_column, 1, 
                                                          ifelse(model_data_spore_filtered$N0_MPN_spore < model_data_spore_filtered$new_column, 0, NA))
model_data_spore_filtered <- model_data_spore_filtered[, !(names(model_data_spore_filtered) %in% c("N0_MPN_spore", "new_column"))]
names(model_data_spore_filtered)[names(model_data_spore_filtered) == "N0_MPN_spore_assigned"] <- "N0_MPN_spore"
filtered_rows <- anti_join(model_data_spore, model_data_spore_filtered, by = c("lot_id", "unit_id"))
model_data_spore <- rbind(filtered_rows,model_data_spore_filtered)
model_data_spore$N0 = model_data_spore$N0_MPN_spore/1900
model_data_spore <- subset(model_data_spore, select = -N0_MPN_spore)
model_data_spore <- as.data.frame(model_data_spore)

# join data
model_data <- rbind(model_data_ppc,model_data_spore)

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
model_data$Nmax <- 10^(model_data$Nmax)

model_data <- model_data %>%
  mutate(Tmin = case_when(
    spoilage_type == "PPC" ~ -4.15,
    spoilage_type == "Spore_Unspoil" ~ -3.62
  ))

model_data$mu_opt <- (model_data$b*(model_data$Topt-model_data$Tmin))^2 

model_data <- as.data.frame(model_data)

# Run simulation  
my_times <- seq(0, 21)
num_iterations <- nrow(model_data)
all_simulations <- list()
for (i in 1:num_iterations) {
  my_primary <- list(mu_opt = model_data$mu_opt[i], Nmax = model_data$Nmax[i], N0 = model_data$N0[i], Q0 = model_data$Q0[i])
  sec_temperature <- list(model = "reducedRatkowsky", xmin = model_data$Tmin[i], b = model_data$b[i], xopt = model_data$Topt[i])
  my_secondary <- list(temperature = sec_temperature)
  growth <- predict_dynamic_growth(times = my_times,
                                   env_conditions = tibble(time = my_times,
                                   temperature = 6),
                                   my_primary,
                                   my_secondary)
  sim <- growth$simulation
  all_simulations[[i]] <- sim
}

final_conc <- do.call(rbind, all_simulations)
df <- final_conc

# Generate output 
model_data_1 <- model_data[rep(1:nrow(model_data), each = 22), ]
model_data_sub <- model_data_1[,c("lot_id","unit_id","N0","spoilage_type", "STorAT")]
count_result <- final_conc[,c("time","N","logN")]
df <- cbind(model_data_sub,count_result)
df$logN[df$logN == -Inf] <- -2

# Validation for PPC
summary_stats_PPC<- list()
for (i in 1:21) {
  filtered_data <- subset(df, time == i)
  percent <- filtered_data %>%
    group_by(lot_id) %>%  
    summarise(percent = sum(logN > 6))
  percent$median_percent <- median(percent$percent)
  percent$perc_5 <- quantile(percent$percent, probs = 0.05)
  percent$perc_95 <- quantile(percent$percent, probs = 0.95)
  percent <- percent[c("median_percent", "perc_5","perc_95")][1,]
  summary_stats_PPC[[i]] <- percent
}
summary_stats_PPC <- do.call(rbind,summary_stats_PPC)
validation_PPC_simulated = summary_stats_PPC[c(7,14,21),]

CVTA_D7_sim = subset(df,time == "7")
summary_stats_CVTA_D7_sim <- quantile(CVTA_D7_sim$logN, probs = c(0.25, 0.5, 0.75))
CVTA_D14_sim = subset(df,time == "14")
summary_stats_CVTA_D14_sim <- quantile(CVTA_D14_sim$logN, probs = c(0.25, 0.5, 0.75))
CVTA_D21_sim = subset(df,time == "21")
summary_stats_CVTA_D21_sim <- quantile(CVTA_D21_sim$logN, probs = c(0.25, 0.5, 0.75))

# Validation for Spore
Microflora <- sample(SPC_D1, 10000, replace = TRUE)
df$logN_correct = log10(df$N + Microflora)

summary_stats_Spore<- list()
for (i in 1:21) {
  filtered_data <- subset(df, time == i)
  percent <- filtered_data %>%
    group_by(lot_id) %>%  
    summarise(percent = sum(logN_correct > 6))
  percent$median_percent <- median(percent$percent)
  percent$perc_5 <- quantile(percent$percent, probs = 0.05)
  percent$perc_95 <- quantile(percent$percent, probs = 0.95)
  percent <- percent[c("median_percent", "perc_5","perc_95")][1,]
  summary_stats_Spore[[i]] <- percent
}
summary_stats_Spore <- do.call(rbind,summary_stats_Spore)
validation_Spore_simulated = summary_stats_Spore[c(7,14,21),]

SPC_D7_sim = subset(df,time == "7")
summary_stats_SPC_D7_sim <- quantile(SPC_D7_sim$logN_correct, probs = c(0.25, 0.5, 0.75))
SPC_D14_sim = subset(df,time == "14")
summary_stats_SPC_D14_sim <- quantile(SPC_D14_sim$logN_correct, probs = c(0.25, 0.5, 0.75))
SPC_D21_sim = subset(df,time == "21")
summary_stats_SPC_D21_sim <- quantile(SPC_D21_sim$logN_correct, probs = c(0.25, 0.5, 0.75))