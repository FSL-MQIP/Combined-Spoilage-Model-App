# Load packages
library(shiny)
library(readxl)
library(tidyverse)
library(purrr)
library(fitdistrplus)
library(EnvStats) 
library(truncnorm)
library(jmuOutlier)
library(ggplot2)

# Load utility functions
source("UtilityFunctions.R")

# Load input files
# AT/ST frequency data
# spore
spore_at_frequency <- read_csv("InputData/spore_ATfreq_clean.csv")
# ppc
ppc_st_frequency <- read_csv("InputData/ppc_STfreq_clean.csv")

# Growth parameter 
# spore
spore_primary_growthparameters <- read_csv("InputData/spore_growthparameters_clean.csv")
# ppc
ppc_primary_growthparameters <- read_csv("InputData/ppc_growthparameters_clean.csv")
# no spoil
noSpoil_growthparameters <- read_csv("InputData/noSpoil_growthparameters.csv")
##join
growth_parameters <-rbind(ppc_primary_growthparameters,spore_primary_growthparameters,noSpoil_growthparameters)

# Set up dataframe for simulation 
n_sim <- 100
n_units <- 10
lot_id <- rep(seq(1, n_sim), each = n_units)
unit_id <- rep(seq(1,n_units), times = n_sim ) 
STorAT <- vector(mode="logical", n_sim  *n_units)

#stage1, storage at facility
t_F <- vector(mode="logical", n_sim  *n_units) 
T_F <- vector(mode="logical", n_sim  *n_units) 
count_F <- vector(mode = "logical", n_sim  *n_units)
#stage2, transport to retail store
t_T <- vector(mode="logical", n_sim  *n_units) 
T_T <- vector(mode="logical", n_sim  *n_units) 
count_T <- vector(mode = "logical", n_sim  *n_units)
#stage3, storage/display at retail store
t_S <- vector(mode="logical", n_sim  *n_units) 
T_S <- vector(mode="logical", n_sim  *n_units) 
count_S <- vector(mode = "logical", n_sim  *n_units)
#stage4, transport from retail store to homes
t_T2 <- vector(mode="logical", n_sim  *n_units) 
T_T2 <- vector(mode="logical", n_sim  *n_units) 
count_T2 <- vector(mode = "logical", n_sim  *n_units)
#stage5, storage at homes
t_H <- vector(mode="logical", n_sim  *n_units) 
T_H <- vector(mode="logical", n_sim  *n_units) 
count_H <- vector(mode = "logical", n_sim  *n_units)
model_data <- data.frame(lot_id, unit_id,STorAT,t_F, T_F, count_F, t_T, T_T, count_T, 
                   t_S, T_S, count_S, t_T2, T_T2, count_T2,t_H, T_H, count_H)

# Define ui
ui <- fluidPage(
  titlePanel("Milk Spoilage Prediction App"),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("Percent_PPC_Spoil","What is the percentage of milk spoiled due to PPC?", value = 40),
      numericInput("count_mean1", "What is the average (mean) PPC bacterial concentration (log10 cfu/mL) in milk?", value = 0.38),
      numericInput("count_sd1", "What is the standard deviation of PPC bacterial concentration (log10 cfu/mL) in milk?", value = 1.11),
      numericInput("Percent_spore_Spoil","What is the percentage of milk spoiled due to sporeformers?", value = 40),
      numericInput("count_mean2", "What is the average (mean) spore concentration (log10 MPN/mL) in milk?", value = -0.72),
      numericInput("count_sd2", "What is the standard deviation of spore concentration (log10 MPN/mL) in milk?", value = 0.99),
      selectInput("threshold","Spoilage threshold", c("US regulation limit (Pasteurized Milk Ordinance): 20,000 CFU/mL" = log10(20000))), 
      numericInput("day","How many storage days at consumer's home to simulate?", value = 35),
      h5(strong("Select spore reduction technology")),
      checkboxInput("mf","Microfiltration (2.2log reduction)", value=FALSE),
      checkboxInput("bf1","Bactofugation single-pass (1.4log reduction)", value=FALSE),
      checkboxInput("bf2","Bactofugation double-pass (2log reduction)", value=FALSE),
      h5(strong("Select PPC reduction technology")),
      checkboxInput("mt1","Improved preventive maintenance (1log reduction)", value=FALSE),
      checkboxInput("mt2","super improved preventive maintenance (3log reduction)", value=FALSE),
      h5(strong("Improved temperature control")),
      checkboxInput("tc","Consumer home storage temperature control (<=4°C)", value=FALSE),
      submitButton("Submit", icon("refresh"))),
    
    mainPanel(
      plotOutput("plot"), 
      verbatimTextOutput("prediction")
      )
    )
  )

# Define server 
server <- function(input, output) {
  model_result <- reactive({
    
# Set seed
set.seed(1)
    
    # Assign temperature profiles 
    # Stage 1: Facility storage
    temps_F <- rep(runif(n_sim,min=3.5,max=4.5),each=n_units) 
    model_data$T_F <- temps_F
    times_F <- rep(runif(n_sim,min=1,max=2),each=n_units)
    model_data$t_F <- times_F
    # Stage 2: Transport from facility to retail
    temps_T <- rep(rtri(n_sim,min=1.7,max=10.0,mode=4.4),each=n_units) 
    model_data$T_T <- temps_T
    times_T <- rep(rtri(n_sim,min=1,max=10,mode=5),each=n_units)
    model_data$t_T <- times_T
    # Stage 3. Retail storage
    unif_mean = 2.3
    unif_b = 5.4
    temps_S <- rep(rtruncnorm(n_sim,a=-1.4,b=unif_b,mean=unif_mean,sd=1.8),each=n_units)
    model_data$T_S <- temps_S
    times_S <- rep(rtruncnorm(n_sim,a=0.042,b=10.0, mean=1.821,sd=3.3),each=n_units)
    model_data$t_S <- times_S
    # Stage 4. Transport from retail to home
    temps_T2 <- rep(rtruncnorm(n_sim,a=0,b=10,mean=8.5,sd=1.0),each=n_units)
    model_data$T_T2 <- temps_T2
    times_T2 <- rep(rtruncnorm(n_sim,a=0.01,b=0.24, mean=0.04,sd=0.02),each=n_units)
    model_data$t_T2 <- times_T2
    # Stage 5. Home storage
    temps <- rep(NA, n_sim)
    for (i in 1:(n_sim*n_units)){
      number <- rlaplace(1,m=4.06,s=2.31)
      while (number > 15 | number < -1) {
        number <- rlaplace(1,m=4.06,s=2.31) #make sure that this cannot be >15 or < -1
      }
      temps[i] <- number
    }
    model_data$T_H <- temps 
    
    # Implement temperature control measures 
    if (input$tc) {
      temp <- rep(NA, n_sim*n_units)
      for (i in 1:(n_sim*n_units)){
        number <- rlaplace(1,m=4.06,s=2.31)
        while (number > 4 | number < -1) {
          number <- rlaplace(1,m=4.06,s=2.31)
        }
        temp[i] <- number
      }
      model_data$T_H <- temps 
    }

# Generate spoilage frequency and assign spoilage types 
num_ppc <- round(input$Percent_PPC_Spoil/ 100 * n_sim*n_units)
num_spore <- round(input$Percent_spore_Spoil/ 100 * n_sim*n_units)
num_no_spoil <- n_sim*n_units - num_ppc - num_spore
model_data$spoilage_type <- c(rep("PPC", num_ppc), rep("Spore", num_spore), rep("No Spoil", num_no_spoil))

# if spoilage type is PPC, multiply time by 24 h (double check)
model_data$t_F <- with(model_data, ifelse(spoilage_type != "PPC", t_F, t_F*24)) 
model_data$t_T <- with(model_data, ifelse(spoilage_type != "PPC", t_T, t_T*24))
model_data$t_S <- with(model_data, ifelse(spoilage_type != "PPC", t_S, t_S*24))
model_data$t_T2 <- with(model_data, ifelse(spoilage_type != "PPC", t_T2, t_T2*24))

# Generate initial count distribution and assign initial count 
# ppc
model_data_ppc <- subset(model_data, spoilage_type == "PPC")
ppc_log10MPN_samp <-  rnorm(n_sim, input$count_mean1, input$count_sd1)
ppc_MPN_samp <- 10^(ppc_log10MPN_samp) 
ppc_MPN_samp_halfgal <- ppc_MPN_samp * 1900
ppc_MPN_init<-vector()
for (i in 1:n_sim){
  ppc_MPN_init_samp <-rep(rpois(n_units, ppc_MPN_samp_halfgal[i]))
  ppc_MPN_init<-c(ppc_MPN_init, ppc_MPN_init_samp)}
ppc_MPN_init[ppc_MPN_init<1] = 0
ppc_MPN_init_mL <- ppc_MPN_init / 1900
ppc_MPN_init_mL <- sample(ppc_MPN_init_mL, num_ppc, replace = T)
model_data_ppc$ppc_MPN_init_mL <- ppc_MPN_init_mL
model_data_ppc$ppc_MPN_init_mL[model_data_ppc$ppc_MPN_init_mL<=0 ] <- 0
model_data_ppc$ppc_MPN_init_mL[model_data_ppc$ppc_MPN_init_mL == 0] <- 0.01
model_data_ppc$ppc_log10_init_mL <- log10(model_data_ppc$ppc_MPN_init_mL)
colnames(model_data_ppc)[ncol(model_data_ppc) - 1] <- "MPN_init_mL"
colnames(model_data_ppc)[ncol(model_data_ppc)] <- "log10_init_mL"

# Add ppc control measures 
if (input$mt1) {model_data_ppc$log10_init_mL = model_data_ppc$log10_init_mL - 1}
if (input$mt2) {model_data_ppc$log10_init_mL = model_data_ppc$log10_init_mL - 3}

# spore
model_data_spore <- subset(model_data, spoilage_type == "Spore")
spore_log10MPN_samp <-  rnorm(n_sim, input$count_mean2, input$count_sd2)
spore_MPN_samp <- 10^(spore_log10MPN_samp) 
spore_MPN_samp_halfgal <- spore_MPN_samp * 1900
spore_MPN_init<-vector()
for (i in 1:n_sim){
  spore_MPN_init_samp <-rep(rpois(n_units, spore_MPN_samp_halfgal[i]))
  spore_MPN_init<-c(spore_MPN_init, spore_MPN_init_samp)}
spore_MPN_init[spore_MPN_init<1] = 0
spore_MPN_init_mL <- spore_MPN_init / 1900
spore_MPN_init_mL <- sample(spore_MPN_init_mL, num_spore, replace = T)
model_data_spore$spore_MPN_init_mL <- spore_MPN_init_mL
model_data_spore$spore_MPN_init_mL[model_data_spore$spore_MPN_init_mL<=0 ] <- 0
model_data_spore$spore_MPN_init_mL[model_data_spore$spore_MPN_init_mL == 0] <- 0.01
model_data_spore$spore_log10_init_mL <- log10(model_data_spore$spore_MPN_init_mL) 
colnames(model_data_spore)[ncol(model_data_spore) - 1] <- "MPN_init_mL"
colnames(model_data_spore)[ncol(model_data_spore)] <- "log10_init_mL"

# Add spore control measures 
if (input$mf) {model_data_spore$log10_init_mL = model_data_spore$log10_init_mL - 2.2}
if (input$bf1) {model_data_spore$log10_init_mL = model_data_spore$log10_init_mL - 1.4}
if (input$bf2) {model_data_spore$log10_init_mL = model_data_spore$log10_init_mL - 2}

# no spoil
df <- read_excel("InputData/initialcounts_060622.xlsx")
noSpoil_initialcounts <- df %>%
  filter(Day_Initial_Actual<=7)%>%
  filter(spoilagetype_actual=="noSpoil")
noSpoil_nfit  <- fitdist(log10(noSpoil_initialcounts$SPC_DI), "norm")
model_data_noSpoil <- subset(model_data, spoilage_type == "No Spoil")
noSpoil_log10MPN_samp <- rnorm(n_sim, noSpoil_nfit$estimate[1], noSpoil_nfit$estimate[2])
noSpoil_MPN_samp <- 10^(noSpoil_log10MPN_samp) 
noSpoil_MPN_samp_halfgal <- noSpoil_MPN_samp * 1900
noSpoil_MPN_init<-vector()
for (i in 1:n_sim){
  noSpoil_MPN_init_samp <-rep(rpois(n_units, noSpoil_MPN_samp_halfgal[i]))
  noSpoil_MPN_init<-c(noSpoil_MPN_init, noSpoil_MPN_init_samp)}
noSpoil_MPN_init[noSpoil_MPN_init<1] = 0
noSpoil_MPN_init_mL <- noSpoil_MPN_init / 1900
noSpoil_MPN_init_mL <- sample(noSpoil_MPN_init_mL, num_no_spoil, replace = T)
model_data_noSpoil$noSpoil_MPN_init_mL <- noSpoil_MPN_init_mL
model_data_noSpoil$noSpoil_MPN_init_mL[model_data_noSpoil$noSpoil_MPN_init_mL<=0 ] <- 0
model_data_noSpoil$noSpoil_MPN_init_mL[model_data_noSpoil$noSpoil_MPN_init_mL == 0] <- 0.01
model_data_noSpoil$noSpoil_log10_init_mL <- log10(model_data_noSpoil$noSpoil_MPN_init_mL)
colnames(model_data_noSpoil)[ncol(model_data_noSpoil) - 1] <- "MPN_init_mL"
colnames(model_data_noSpoil)[ncol(model_data_noSpoil)] <- "log10_init_mL"

# Assign AT/ST types 
# ppc
model_data_ppc$STorAT <- sample(ppc_st_frequency$ClosestST, num_ppc,replace=T)
model_data_ppc$STorAT <- paste0("ST_", model_data_ppc$STorAT)
# spore
model_data_spore$STorAT <- sample(spore_at_frequency$ClosestAT, num_spore,replace=T)
model_data_spore$STorAT <- paste0("AT_", model_data_spore$STorAT)
# no spoil
model_data_noSpoil$STorAT <- rep("NS_0", num_no_spoil)
# join data
model_data <- rbind(model_data_ppc,model_data_spore,model_data_noSpoil)

# Generate allele index
model_data$allele_index <- match(model_data$STorAT, growth_parameters$STorAT)

# Add log10Nmax
model_data$log10Nmax = growth_parameters$LOG10Nmax[model_data$allele_index]

# Add model type 
model_data$model_type = growth_parameters$model_type[model_data$allele_index]

# Simulation 
# Stage 1: Facility storage 
# Determine NewLag_F and NewMu_F
model_data$newLag_F <- lagAtNewTemp(model_data$T_F, growth_parameters$lag[model_data$allele_index], T0=growth_parameters$T0[model_data$allele_index])
model_data$newMu_F <- muAtNewTemp(model_data$T_F,growth_parameters$mumax[model_data$allele_index],T0=growth_parameters$T0[model_data$allele_index])
# Determine count_F
model_data = model_data %>%
  rowwise() %>% 
  mutate(count_F = log10N(log10_init_mL, log10Nmax, newLag_F, newMu_F, t_F, model_type))

# Stage 2: Transport from facility to retail store
# Determine NewLag_T and NewMu_T
model_data$newLag_T <- lagAtNewTemp(model_data$T_T, growth_parameters$lag[model_data$allele_index], T0=growth_parameters$T0[model_data$allele_index])
model_data$newMu_T <- muAtNewTemp(model_data$T_T,growth_parameters$mumax[model_data$allele_index],T0=growth_parameters$T0[model_data$allele_index])
# Determine count_T
model_data = model_data %>% 
  rowwise() %>% 
  mutate(Lag_T = max(0, 1 - t_F/newLag_F) * newLag_T) %>%  
  mutate(Mu_T = if_else(T_T >= T_F * 0.75 & T_T <= T_F * 1.25, newMu_F, newMu_T)) %>% 
  mutate(count_T = log10N(count_F, log10Nmax, Lag_T, Mu_T, t_T, model_type))

# Stage 3: Display at retail
# Determine NewLag_S and NewMu_S
model_data$newLag_S <- lagAtNewTemp(model_data$T_S, growth_parameters$lag[model_data$allele_index], T0=growth_parameters$T0[model_data$allele_index])
model_data$newMu_S <- muAtNewTemp(model_data$T_S,growth_parameters$mumax[model_data$allele_index],T0=growth_parameters$T0[model_data$allele_index])
# Determine count_S
model_data = model_data %>% 
  rowwise() %>% 
  mutate(Lag_S = max(0, 1 - t_T/Lag_T) * newLag_S) %>%  
  mutate(Mu_S = if_else(T_S >= T_T * 0.75 & T_S <= T_T * 1.25, newMu_T, newMu_S)) %>% 
  mutate(count_S = log10N(count_T, log10Nmax, Lag_S, Mu_S, t_S, model_type))

# Stage 4: Transport from retail store to consumer home
# Determine NewLag_T2 and NewMu_T2
model_data$newLag_T2 <- lagAtNewTemp(model_data$T_T2, growth_parameters$lag[model_data$allele_index], T0=growth_parameters$T0[model_data$allele_index])
model_data$newMu_T2 <- muAtNewTemp(model_data$T_T2,growth_parameters$mumax[model_data$allele_index],T0=growth_parameters$T0[model_data$allele_index])
# Determine count_T2
model_data = model_data %>% 
  rowwise() %>% 
  mutate(Lag_T2 = max(0, 1 - t_S/Lag_S) * newLag_T2) %>%  
  mutate(Mu_T2 = if_else(T_T2 >= T_S * 0.75 & T_T2 <= T_S * 1.25, newMu_S, newMu_T2)) %>% 
  mutate(count_T2 = log10N(count_S, log10Nmax, Lag_T2, Mu_T2, t_T2, model_type))

# Stage 5: Storage at homes
# Determine NewLag_H and NewMu_H
model_data$newLag_H <- lagAtNewTemp(model_data$T_H, growth_parameters$lag[model_data$allele_index], T0=growth_parameters$T0[model_data$allele_index])
model_data$newMu_H <- muAtNewTemp(model_data$T_H,growth_parameters$mumax[model_data$allele_index],T0=growth_parameters$T0[model_data$allele_index])
# Determine count_H
end_day = input$day
model_data = model_data %>% 
  slice(rep(1:n(), each = end_day)) 

model_data$t_H = rep(1:end_day, times = n_sim * n_units)
model_data$t_H <- with(model_data, ifelse(spoilage_type != "PPC", t_H, t_H*24))

model_data = model_data %>% 
  rowwise() %>% 
  mutate(Lag_H = max(0, 1 - t_T2/Lag_T2) * newLag_H) %>%  
  mutate(Mu_H = if_else(T_H >= T_T2 * 0.75 & T_H <= T_T2 * 1.25, newMu_T2, newMu_H)) %>% 
  mutate(count_H = log10N(count_T2, log10Nmax, Lag_H, Mu_H, t_H, model_type))

# Determine percent spoiled
model_data$t_H = rep(1:end_day, times = n_sim * n_units)
filtered_counts <- map(1:end_day, ~ filter(model_data, t_H == .x) %>% pull(count_H))
percentages <- map(filtered_counts, ~ sum(. > input$threshold) / length(.))
percent_spoil <- data.frame(day = 1:end_day, percentage = unlist(percentages))
percent_spoil$percentage = percent_spoil$percentage * 100
percent_spoil
})

# Generate Output
  # Plot
  output$plot <- renderPlot({
    ggplot(data = model_result(), aes(x = day, y = percentage))+
      geom_line(aes(y=percentage))+
      labs(title="Simulated % of spoiled half-gallon milk containers",
           x="Consumer storage (days)",
           y="% of spoiled half-gallon milk containers")+
      theme_classic()+
      ylim(0,100) + 
      scale_x_continuous(breaks = scales::pretty_breaks(10))
  })
  
  # Text
  output$prediction <- renderText({
    df = model_result()
    subset_df <- subset(df, percentage > 50)
    shelflife <- subset_df$day[1]
    text <- paste("The predicted shelf life is", shelflife, "days. The shelf life is defined as the first day when over 50% milk containers exceed the spoilage threshold.")
    text
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

