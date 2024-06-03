training_data <- read.csv("train_data_with_DI.csv")
testing_data <- read.csv("test_data_with_DI_new.csv")
combined_data <- rbind(training_data[,3:14],testing_data[,c(2:5,7,9,11,6,8,10,12:13)])
combined_data_sub <- subset(combined_data, Day_Initial <= 3)

combined_data_PPC <- subset(combined_data_sub,spoilagetype == "PPC")
CVTA_D7 <- log10(combined_data_PPC$CVTA_D7)
CVTA_D7 <- na.omit(CVTA_D7)
CVTA_D14 <- log10(combined_data_PPC$CVTA_D14)
CVTA_D14 <- na.omit(CVTA_D14)
CVTA_D21 <- log10(combined_data_PPC$CVTA_D21)
CVTA_D21 <- na.omit(CVTA_D21)

combined_data_Spore <- subset(combined_data_sub,spoilagetype %in% c("spore spoilage", "no spoilage"))
SPC_D1 <- combined_data_sub$SPC_DI
SPC_D7 <- log10(combined_data_Spore$SPC_D7)
SPC_D7 <- na.omit(SPC_D7)
SPC_D14 <- log10(combined_data_Spore$SPC_D14)
SPC_D14 <- na.omit(SPC_D14)
SPC_D21 <- log10(combined_data_Spore$SPC_D21)
SPC_D21 <- na.omit(SPC_D21)

# Validate PPC using Kolmogorov-Smirnov test
ks_result_D7_ppc <- ks.test(CVTA_D7_sim, CVTA_D7)
ks_result_D14_ppc <- ks.test(CVTA_D14_sim, CVTA_D14)
ks_result_D21_ppc <- ks.test(CVTA_D21_sim, CVTA_D21)

# Validate Sporeformers using Kolmogorov-Smirnov test
ks_result_D7_spore <- ks.test(SPC_D7_sim, SPC_D7)
ks_result_D14_spore <- ks.test(SPC_D14_sim, SPC_D14)
ks_result_D21_spore <- ks.test(SPC_D21_sim, SPC_D21)
