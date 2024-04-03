    # Set work directory 
    setwd("C:/Users/sujun/Documents/GitHub/Combined-Spoilage-Model-App")
    
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
    no_spoil_initial <- read.csv("InputData/initialcounts_060622.csv")
    # Growth parameters
    data_ppc <- read.csv("InputData/ppc_gp.csv")
    data_ppc$mumax <- 0.684 * data_ppc$mumax
    data_sporeformer <- read.csv("InputData/sporeformer_gp.csv")
    gp_noSpoil <- read.csv("InputData/NoSpoil_gp.csv")
    colnames(gp_noSpoil)[1] <- "isolate"
    # AT frequency
    ppc_ST_freq <- read.csv("InputData/ppc_STfreq_clean.csv")
    spore_AT_freq <- read.csv("InputData/spore_ATfreq_new.csv")
    
    # Set seed
    set.seed(1)
    
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
    # data$T_F <- rep(runif(n_sim*n_unit,min=3.5,max=4.5)) #uniform distribution
    ## (b) Sample the storage time (in days) distribution
    # data$t_F <- rep(runif(n_sim*n_unit,min=1,max=2)) #uniform distribution
    
    # Stage 2: transport from facility to retail store
    ## (a)  Sample the temperature distribution
    # data$T_T <- rep(rtri(n_sim*n_unit,min=1.7,max=10.0,mode=4.4)) #triangular distribution
    ## (b) Sample the transportation time (in days) distribution
    # data$t_T <- rep(rtri(n_sim*n_unit,min=1,max=10,mode=5))
    
    # Stage 3: storage/display at retail store
    ## (a)  Sample the temperature distribution
    # data$T_S <- rep(rtruncnorm(n_sim*n_unit,a=-1.4,b=5.4,mean=2.3,sd=1.8)) #truncated normal distribution
    ## (b) Sample the storage time (in days) distribution
    # data$t_S <- rep(rtruncnorm(n_sim*n_unit,a=0.042,b=10.0, mean=1.821,sd=3.3)) #truncated normal distribution
    
    ## Stage 4: transportation from retail store to home
    ## (a)  Sample the temperature distribution
    # data$T_T2 <- rep(rtruncnorm(n_sim*n_unit,a=0,b=10,mean=8.5,sd=1.0)) #truncated normal distribution
    ## (b) Sample the transportation time (in days) distribution 
    # data$t_T2 <- rep(rtruncnorm(n_sim*n_unit,a=0.01,b=0.24, mean=0.04,sd=0.02)) #truncated normal distribution
    
    ## Stage 5: home storage 
    ## (a)  Sample the temperature distribution
    # temps <- rep(NA, n_sim*n_unit)
    # for (i in 1:(n_sim*n_unit)){
      # number <- rlaplace(1,m=4.06,s=2.31)
      # while (number > 15 | number < -1) {
        # number <- rlaplace(1,m=4.06,s=2.31) #truncated laplace distribution 
      # }
      # temps[i] <- number
    # }
    # data$T_H <- temps
    
    ## (b) Define shelf-life day for all units
    ## Day 35
    # data$t_H = 35
    
    # Generate spoilage frequency and assign spoilage types
    # PPC Spoilage % (PPC+)  
    # Good Plant
    # data <- data %>%
      # group_by(lot_id) %>%
      # mutate (P_ppc = runif(1, 0.125, 0.313))
    
    # Medium Plant
    # data <- data %>%
      # group_by(lot_id) %>%
      # mutate(P_ppc = runif(1, 0.367, 0.667))
    
    # Bad Plant
    # data <- data %>%
      # group_by(lot_id) %>%
      # mutate (P_ppc = runif(1, 0.75, 1))
    
    # Spore Spoilage % (PPC- & Spore+)
    # data <- data %>%
      # group_by(lot_id) %>%
      # mutate(P_spore = (1 - P_ppc) * runif(1, 0.37, 0.502))
    
    # Sanity check 
    data$P_ppc = 0
    data$P_spore = 1
    
    # Assign spoilage type
    result_list <- list()
    for (i in 1:n_sim) {
      data_lot <- subset(data, lot_id == as.character(i))
      num_ppc <- round(unique(data_lot$P_ppc) * n_unit)
      num_spore <- round(unique(data_lot$P_spore) * n_unit)
      num_no_spoil <- n_unit - num_ppc - num_spore
      spoiler_types <- c(rep("PPC", num_ppc), rep("Spore", num_spore), rep("No Spoil", num_no_spoil))
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
    model_data_spore <- subset(data, spoilage_type == "Spore")
    model_data_spore$STorAT <- rep(15, nrow(model_data_spore))
    # model_data_spore$STorAT <- sample(spore_AT_freq$ClosestAT, nrow(model_data_spore),replace=T)
    model_data_spore$STorAT <- paste0("AT_", model_data_spore$STorAT)
    # no spoil
    model_data_noSpoil <- subset(data, spoilage_type == "No Spoil")
    model_data_noSpoil$STorAT <- rep("NS_0", nrow(model_data_noSpoil))
    
    # Assign initial contamination distributions
    # ppc
    model_data_ppc <- model_data_ppc %>%
      rowwise() %>%
      mutate(N0 = rnorm(n = 1, mean = 0.38, sd = 1.11))
    model_data_ppc <- as.data.frame(model_data_ppc)
    
    # spore 
    model_data_spore <- model_data_spore %>%
      rowwise() %>%
      mutate(N0 = rnorm(n = 1, mean = -0.72, sd = 0.99))
    model_data_spore <- as.data.frame(model_data_spore)
    
    # no Spoil
    no_spoil_initial <- no_spoil_initial %>%
      filter(Day_Initial_Actual<=7)%>%
      filter(spoilagetype_actual=="noSpoil")
    noSpoil_nfit  <- fitdist(log10(no_spoil_initial$SPC_DI), "norm")
    initial_mean <- noSpoil_nfit$estimate[1]
    initial_sd <- noSpoil_nfit$estimate[2]
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
    model_data$Nmax <- 10^(model_data$Nmax)
    
    model_data$N0 <- 10^(model_data$N0)
    
    model_data <- model_data %>%
      mutate(Tmin = case_when(
        spoilage_type == "PPC" ~ -4.15,
        spoilage_type == "Spore" ~ -3.62,
        spoilage_type == "No Spoil" ~ 0
      ))
    
    model_data$mu_opt <- (model_data$b*(model_data$Topt-model_data$Tmin))^2 
    
    model_data <- model_data %>%
      mutate(mu_opt = ifelse(spoilage_type == "No Spoil", 0, mu_opt),
             Q0 = ifelse(spoilage_type  == "No Spoil", 0, Q0),
             b = ifelse(spoilage_type == "No Spoil", 0, b))
    
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
    my_times <- seq(0, 21, length = 22)
    num_iterations <- nrow(model_data)
    all_simulations <- list()
    for (i in 1:num_iterations) {
      my_primary <- list(mu_opt = model_data$mu_opt[i], Nmax = model_data$Nmax[i], N0 = model_data$N0[i], Q0 = model_data$Q0[i])
      sec_temperature <- list(model = "reducedRatkowsky", xmin = model_data$Tmin[i], b = model_data$b[i], xopt = model_data$Topt[i])
      my_secondary <- list(temperature = sec_temperature)
      growth <- predict_dynamic_growth(times = my_times,
                                       # env_conditions = tibble(time = env_cond_time[i,],
                                                               # temperature = env_cond_temp[i,]),
                                       env_conditions = tibble(time = my_times,
                                                               temperature = 6),
                                       my_primary,
                                       my_secondary)
      sim <- growth$simulation
      all_simulations[[i]] <- sim
     }
    
    final_conc <- do.call(rbind, all_simulations)
    final_conc$logN[is.na(final_conc$logN)] <- model_data_noSpoil$N0
    final_conc$N[is.na(final_conc$N)] <- 10^(model_data_noSpoil$N0)
    
    # Generate output 
    model_data_1 <- model_data[rep(1:nrow(model_data), each = 22), ]
    model_data_sub <- model_data_1[,c("lot_id","unit_id","N0","spoilage_type", "STorAT")]
    # model_data_sub <- model_data_1[,c("lot_id","unit_id","spoilage_type", "t_F", "T_F", "t_T", "T_T", 
                                    # "t_S", "T_S", "t_T2", "T_T2", "T_H")]
    count_result <- final_conc[,c("time","N","logN")]
    df <- cbind(model_data_sub,count_result)
    
    df_D7 = subset(df, time == "7")
    df_D14 = subset(df, time == "14")
    df_D21 = subset(df, time == "21")
    
    # microflora
    Microflora <- sample(SPC_D1, 10000, replace = TRUE)
    
    Spore_D7 =  df_D7$N 
    SPC_D7_sim = log10(Spore_D7 + Microflora)
    # CVTA_D7_sim = df_D7$logN
    
    Spore_D14 = df_D14$N
    SPC_D14_sim = log10(Spore_D14 + Microflora)
    # CVTA_D14_sim = df_D14$logN
    
    Spore_D21 = df_D21$N
    SPC_D21_sim = log10(Spore_D21 + Microflora)
    # CVTA_D21_sim = df_D21$logN
    
    # percent_spoiled_spore_6dC <- df %>%
      # group_by(time) %>%
      # summarise(percent = sum(logN > log10(20000))/10000*100)
    
    # Calculate the mean and median conc. for each day 
    result_bad <- df %>%
      group_by(time) %>%
      summarize(mean_logN = mean(logN),
                median_logN = median(logN),
                meanN = mean(N),
                log_meanN = log10(mean(N)))
    
    # Calculate the sum of logN over 6 logs for each lot_id
    summary_stats_bad <- list()
    for (i in 1:35) {
      filtered_data <- subset(df, time == i)
      percent <- filtered_data %>%
        group_by(lot_id) %>%
        summarise(percent = sum(logN > 6))
      percent$mean_percent <- mean(percent$percent)
      percent$median_percent <- median(percent$percent)
      percent$perc_2.5 <- quantile(percent$percent, probs = 0.025)
      percent$perc_97.5 <- quantile(percent$percent, probs = 0.975)
      percent <- percent[c("mean_percent", "median_percent", "perc_2.5", "perc_97.5")][1,]
      summary_stats_bad[[i]] <- percent
    }
    summary_stats_bad <- do.call(rbind,summary_stats_bad)
    
write.csv(summary_stats_good, "summary_stats_good_new.csv")
write.csv(result_good, "result_good_new.csv")
write.csv(summary_stats_medium, "summary_stats_medium_new.csv")
write.csv(result_medium, "result_medium_new.csv")
write.csv(summary_stats_bad, "summary_stats_bad_new.csv")
write.csv(result_bad, "result_bad_new.csv")

# plot 
ggplot(spore_count_6dC, aes(x = logN, fill = factor(time))) +
  geom_histogram(binwidth = 0.5, position = "dodge", color = "black") +
  facet_wrap(~ time, ncol = 3) +
  geom_vline(xintercept = 4.3, linetype = "dashed", color = "red", size = 1) +
  scale_x_continuous(breaks = seq(floor(min(spore_count_6dC$logN)), ceiling(max(spore_count_6dC$logN)), by = 1)) +
  labs(title = "Spore_Count_6dC",
       x = "logN",
       y = "Frequency") +
  theme_minimal()
