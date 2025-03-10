    # Set work directory 
    setwd("C:/Users/sujun/Documents/GitHub/Combined-Spoilage-Model-App")
    
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
    temps <- rep(NA, n_sim*n_unit)
    for (i in 1:(n_sim*n_unit)){
      number <- rlaplace(1,m=4.06,s=2.31)
      while (number > 15 | number < -1) {
        number <- rlaplace(1,m=4.06,s=2.31) #truncated laplace distribution 
      }
      temps[i] <- number
    }
    data$T_H <- temps
    
    # scenario 4) consumer home temperature <= 4C 
    # temps <- rep(NA, n_sim*n_unit)
    # for (i in 1:(n_sim*n_unit)){
    # number <- rlaplace(1,m=4.06,s=2.31)
    # while (number > 4 | number < -1) {
    # number <- rlaplace(1,m=4.06,s=2.31)
    # }
    # temps[i] <- number
    # }
    # data$T_H <- temps
    
    ## (b) Define shelf-life day for all units
    data$t_H = 28
    
    # scenario 5) and 6) transportation temperature between facility and retail <= 7C or 5C (run this the last)
    # data$T_T <- ifelse(data$T_T > 7, {
      # repeat {
        # new_value <- rtri(1, min = 1.7, max = 10.0, mode = 4.4)
        # if (new_value <= 5) break
      # }
      # new_value
    # }, data$T_T)
    
    # scenario 7) transportation time between facility and retail <= 2 days
    data$t_T <- ifelse(data$t_T > 3, {
      repeat {
        new_value <- rtri(1, min = 1, max = 10, mode = 5)
        if (new_value <= 3) break
      }
      new_value
    }, data$t_T)
    
    # PPC Spoilage %
    # Good Plant
    # data <- data %>%
      # group_by(lot_id) %>%
      # mutate (P_ppc = runif(1, 0.125, 0.294))
    
    # Medium Plant
    # data <- data %>%
      # group_by(lot_id) %>%
      # mutate(P_ppc = runif(1, 0.4, 0.625))
    
    # Bad Plant
    data <- data %>%
      group_by(lot_id) %>%
      mutate (P_ppc = runif(1, 0.778, 1))
    
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
    
    # scenario 3) reduce initial ppc by 1-3 log
    # log_reduction_ppc = rep(runif(n_sim,min=1,max=3))
    # log_reduction_ppc = data.frame(lot_id = c(1:100), log_reduction = log_reduction_ppc)
    # model_data_ppc <- model_data_ppc %>%
      # left_join(log_reduction_ppc, by = "lot_id")
    # model_data_ppc$N0 = model_data_ppc$N0 - model_data_ppc$log_reduction
    # model_data_ppc <- subset(model_data_ppc, select = -log_reduction)
    
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
    
    # scenario 1) and 2) reduce initial spore by 0.6-3.1 log or 1.23-1.64 log
    # log_reduction_spore = rep(runif(n_sim,min=1.23,max=1.64))
    # log_reduction_spore = data.frame(lot_id = c(1:100), log_reduction = log_reduction_spore)
    # model_data_spore <- model_data_spore %>%
    # left_join(log_reduction_spore, by = "lot_id")
    # model_data_spore$N0 = model_data_spore$N0 - model_data_spore$log_reduction
    # model_data_spore <- subset(model_data_spore, select = -log_reduction)
    
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
    
    # Model temperature profiles of 10000 units HTST milk
    env_cond_time <- matrix(c(rep(0,10000),
                              model_data$t_F, 
                              model_data$t_F+0.00001,
                              model_data$t_F + model_data$t_T,
                              model_data$t_F + model_data$t_T+0.00001,
                              model_data$t_F + model_data$t_T + model_data$t_S,
                              model_data$t_F + model_data$t_T + model_data$t_S+0.00001,
                              model_data$t_F + model_data$t_T + model_data$t_S + model_data$t_T2,
                              model_data$t_F + model_data$t_T + model_data$t_S + model_data$t_T2+0.00001,
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
    my_times <- seq(0, 28)
    num_iterations <- nrow(model_data)
    all_simulations <- list()
    for (i in 1:num_iterations) {
      my_primary <- list(mu_opt = model_data$mu_opt[i], Nmax = model_data$Nmax[i], N0 = model_data$N0[i], Q0 = model_data$Q0[i])
      sec_temperature <- list(model = "reducedRatkowsky", xmin = model_data$Tmin[i], b = model_data$b[i], xopt = model_data$Topt[i])
      my_secondary <- list(temperature = sec_temperature)
      growth <- predict_dynamic_growth(times = my_times,
                                       env_conditions = tibble(time = env_cond_time[i,],
                                                               temperature = env_cond_temp[i,]),
                                       # env_conditions = tibble(time = my_times,
                                                               # temperature = 6),
                                       my_primary,
                                       my_secondary)
      sim <- growth$simulation
      all_simulations[[i]] <- sim
     }
    
    final_conc <- do.call(rbind, all_simulations)
    df <- final_conc
    
    # Generate output 
    model_data_1 <- model_data[rep(1:nrow(model_data), each = 29), ]
    model_data_sub <- model_data_1[,c("lot_id","unit_id","N0","spoilage_type", "STorAT")]
    count_result <- final_conc[,c("time","N","logN")]
    df <- cbind(model_data_sub,count_result)
    df$logN[df$logN == -Inf] <- -2
    
    # Get output for sensitivity analysis
    # model_data_sub <- model_data_1[,c("t_F", "T_F", "t_T", "T_T", "t_S", "T_S", 
                                      # "t_T2", "T_T2", "T_H","N0")]
    # count_result <- final_conc[,c("time","logN")]
    # df <- cbind(model_data_sub,count_result)
    # df$logN[df$logN == -Inf] <- -2
    # Good_Plant_D21_SA = subset(df, time == "21")
    # Medium_Plant_D21_SA = subset(df, time == "21")
    # Medium_Plant_D14_SA = subset(df, time == "14")
    # Bad_Plant_D21_SA = subset(df, time == "21")
    # Bad_Plant_D7_SA = subset(df, time == "7")
    # Good_Plant_D21_SA$time <- NULL
    # Medium_Plant_D21_SA$time <- NULL
    # Medium_Plant_D14_SA$time <- NULL
    # Bad_Plant_D21_SA$time <- NULL
    # Bad_Plant_D7_SA$time <- NULL
    # write.csv(Good_Plant_D21_SA, "Good_Plant_D21_SA.csv")
    # write.csv(Medium_Plant_D21_SA, "Medium_Plant_D21_SA.csv")
    # write.csv(Medium_Plant_D14_SA, "Medium_Plant_D14_SA.csv")
    # write.csv(Bad_Plant_D21_SA, "Bad_Plant_D21_SA.csv")
    # write.csv(Bad_Plant_D7_SA, "Bad_Plant_D7_SA.csv")
    
    # Baseline model output (percent > 6 log over shelf-life)
    # Calculate the sum of logN over 6 logs for each lot_id
    summary_stats_Good_Total <- list()
    for (i in 1:28) {
      filtered_data <- subset(df, time == i)
      percent <- filtered_data %>%
        group_by(lot_id) %>%  
        summarise(percent = sum(logN > 6))
      percent$median_percent <- median(percent$percent)
      percent$perc_5 <- quantile(percent$percent, probs = 0.05)
      percent$perc_95 <- quantile(percent$percent, probs = 0.95)
      percent <- percent[c("median_percent", "perc_5","perc_95")][1,]
      summary_stats_Good_Total[[i]] <- percent
    }
    summary_stats_Good_Total <- do.call(rbind,summary_stats_Good_Total)
    summary_stats_Good_Total$Plant = "Good"
    summary_stats_Good_Total$Bacteria = "Total"
    summary_stats = rbind(summary_stats_Good_Total, summary_stats_Good_PPC, summary_stats_Good_Spore, 
                          summary_stats_Medium_Total, summary_stats_Medium_PPC, summary_stats_Medium_Spore, 
                          summary_stats_Bad_Total, summary_stats_Bad_PPC, summary_stats_Bad_Spore)
    summary_stats$day <- rep(1:28, times = 9)
    
    # Figure 2
    summary_stats$Plant <- factor(summary_stats$Plant, levels = c("Good", "Medium", "Bad"))
    
    custom_labels <- c("Good" = "Long Shelf-life Milk Plant",
                       "Medium" = "Medium Shelf-life Milk Plant",
                       "Bad" = "Short Shelf-life Milk Plant")
    
    bacteria_labels <- c("Total" = "Total Bacterial Count",
                         "PPC" = "Gram-negative Bacteria",
                         "Spore" = "Spore-forming Bacteria")
    
    perc_spoil = ggplot(summary_stats, aes(x = day)) +
      geom_line(aes(y = median_percent, color = Bacteria, linetype = "Median"), size = 1) +
      geom_line(aes(y = perc_5, color = Bacteria, linetype = "5th and 95th Percentiles"), size = 0.5) +
      geom_line(aes(y = perc_95, color = Bacteria, linetype = "5th and 95th Percentiles"), size = 0.5) +
      geom_hline(yintercept = 25, linetype = "solid", color = "red") +
      
      facet_wrap(~ Plant, labeller = as_labeller(custom_labels), scales = "free_y") +  # Allow free y-scale and space facets more
      labs(x = "Days of refrigerated storage", 
           y = "Proportion of spoiled milk containers",
           color = "Bacteria Type") +
      guides(linetype = "none") +
      scale_x_continuous(limits = c(0, 28), breaks = seq(0, 28, 4), expand = c(0, 0)) +  # Keep no padding around x-axis
      scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), labels = scales::percent_format(scale = 1), expand = c(0, 0)) +  # Keep no padding around y-axis
      scale_color_manual(values = c("Total" = "blue", 
                                    "PPC" = "brown", 
                                    "Spore" = "darkgreen"),
                         labels = bacteria_labels) +
      scale_linetype_manual(values = c("Median" = "solid", "5th and 95th Percentiles" = "dashed")) +
      theme_minimal() +
      theme(panel.grid.major = element_blank(),  
            panel.grid.minor = element_blank(),
            axis.ticks = element_line(color = "black"), 
            axis.ticks.length = unit(0.2, "cm"),
            axis.line = element_line(color = "black"),
            strip.background = element_blank(),  
            strip.text = element_text(size = 14, margin = margin(b = 10)),
            axis.title = element_text(size = 14),  
            axis.text = element_text(size = 12),  
            legend.text = element_text(size = 12))  
    
    ggsave(filename = "Figure2_final.png", 
           plot = perc_spoil, 
           width = 15, height = 7, units = "in", dpi = 300,
           bg = "white")
    
    # Figure 3
    df_Bad_D21 <- subset(df, time == "21")
    percentiles_Bad_D21 <- df_Bad_D21 %>%
      group_by(lot_id) %>%
      summarize(
        p0 = quantile(logN, 0),
        p10 = quantile(logN, 0.10),
        p20 = quantile(logN, 0.20),
        p30 = quantile(logN, 0.30),
        p40 = quantile(logN, 0.40),
        p50 = quantile(logN, 0.50),
        p60 = quantile(logN, 0.60),
        p70 = quantile(logN, 0.70),
        p80 = quantile(logN, 0.80),
        p90 = quantile(logN, 0.90),
        p100 = quantile(logN, 1)
      )
    medians_Bad_D21 <- apply(percentiles_Bad_D21[-1], 2, median)
    percentile_5_Bad_D21 <- apply(percentiles_Bad_D21[-1], 2, function(x) quantile(x, 0.05))
    percentile_95_Bad_D21 <- apply(percentiles_Bad_D21[-1], 2, function(x) quantile(x, 0.95))
    result_Bad_D21 = rbind(medians_Bad_D21, percentile_5_Bad_D21, percentile_95_Bad_D21)
    quantiles = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)
    result_Bad_D21 = rbind(quantiles,result_Bad_D21)
    rownames(result_Bad_D21) = c("quantiles","medians","percentile_5","percentile_95")
    result_Bad_D21 = t(result_Bad_D21)
    result_Bad_D21 = as.data.frame(result_Bad_D21)
    
    datasets <- list(result_Good_D7, result_Medium_D7,result_Bad_D7,
                     result_Good_D14, result_Medium_D14, result_Bad_D14, 
                     result_Good_D21,result_Medium_D21, result_Bad_D21)
    
    titles <- c("(A) Long Shelf-life Milk Plant Day 7", "(D) Medium Shelf-life Milk Plant Day 7", "(G) Short Shelf-life Milk Plant Day 7",
                "(B) Long Shelf-life Milk Plant Day 14", "(E) Medium Shelf-life Milk Plant Day 14", "(H) Short Shelf-life Milk Plant Day 14",
                "(C) Long Shelf-life Milk Plant Day 21", "(F) Medium Shelf-life Milk Plant Day 21", "(I) Short Shelf-life Milk Plant Day 21")
    
    plots <- lapply(1:length(datasets), function(i) {
      ggplot(datasets[[i]], aes(x = medians, y = quantiles)) +
        geom_line(color = "black") +
        geom_line(aes(x = percentile_5, y = quantiles), color = "black", linetype = "dashed") +
        geom_line(aes(x = percentile_95, y = quantiles), color = "black", linetype = "dashed") +
        geom_vline(xintercept = 6, color = "black", linetype = "solid") + 
        xlim(-2, 8) +
        ylim(0, 1) +
        labs(x = expression("log"[10]*" CFU/mL"), y = "Percentiles", title = titles[i]) +
        scale_x_continuous(breaks = seq(-2, 8, by = 1)) +  
        scale_y_continuous(breaks = seq(0, 1, by = 0.1)) + 
        theme_minimal() +
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.ticks = element_line(color = "black", size = 0.5), 
              axis.line = element_line(color = "black", size = 0.5))  
    })
    bacterial_concentration =  
      grid.arrange(plots[[1]], plots[[2]], plots[[3]], 
                   plots[[4]], plots[[5]], plots[[6]], 
                   plots[[7]], plots[[8]], plots[[9]], 
                   nrow = 3, ncol = 3)
    ggsave(filename = "Figure3_final.png", 
           plot = bacterial_concentration, 
           width = 15, height = 7, units = "in", dpi = 300,
           bg = "white")
    
