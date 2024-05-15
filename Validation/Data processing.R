combined_data <- read.csv("combined_data.csv")
combined_data$SampleID <- gsub("-\\d+", "", combined_data$SampleID)
filtered_combined_data <- subset(combined_data, !(endsWith(SampleID, "5") | endsWith(SampleID, "6")))
write.csv(filtered_combined_data,"filtered_combined_data.csv")

combined_data_sub <- subset(filtered_combined_data, Day_Initial <= 3)

combined_data_PPC <- subset(combined_data_sub,spoilagetype == "PPC")
CVTA_D7 <- log10(combined_data_PPC$CVTA_D7)
CVTA_D7 <- na.omit(CVTA_D7)
CVTA_D14 <- log10(combined_data_PPC$CVTA_D14)
CVTA_D14 <- na.omit(CVTA_D14)
CVTA_D21 <- log10(combined_data_PPC$CVTA_D21)
CVTA_D21 <- na.omit(CVTA_D21)

combined_data_Spore <- subset(combined_data_sub, spoilagetype %in% c("spore spoilage", "no spoilage"))
SPC_D1 <- combined_data_sub$SPC_DI
SPC_D7 <- log10(combined_data_Spore$SPC_D7)
SPC_D7 <- na.omit(SPC_D7)
SPC_D14 <- log10(combined_data_Spore$SPC_D14)
SPC_D14 <- na.omit(SPC_D14)
SPC_D21 <- log10(combined_data_Spore$SPC_D21)
SPC_D21 <- na.omit(SPC_D21)

hist(SPC_D7, main = "Histogram of SPC Count (VSL D7)", xlab = "log10CFU/mL", breaks = 25)
hist(SPC_D14, main = "Histogram of SPC Count (VSL D14)", xlab = "log10CFU/mL", breaks = 25)
hist(SPC_D21, main = "Histogram of SPC Count (VSL D21)", xlab = "log10CFU/mL", breaks = 25)

