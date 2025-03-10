library(dplyr)
library(tidyr)
library(ggplot2)

test_data = read.csv("test_data_tier_plant.csv")
train_data = read.csv("train_data_tier_plant.csv")
combined_dataset = rbind(test_data,train_data)
combined_dataset$SampleID_real <- gsub("-\\d+", "", combined_dataset$SampleID_real)
# Filter out flavor milk samples 
filtered_combined_data <- subset(combined_dataset, !(endsWith(SampleID_real, "5") | endsWith(SampleID_real, "6")))
filtered_combined_data <- filtered_combined_data[,c(5,6,17)]

# Not group by Tier 1, 2, 3
# Group by plantID
combined_dataset_grouped <- filtered_combined_data %>%
  group_by(plantID, spoilagetype) %>%
  summarise(count = n()) %>%
  ungroup()
combined_dataset_grouped <- combined_dataset_grouped %>%
  pivot_wider(names_from = spoilagetype, values_from = count, values_fill = 0)
combined_dataset_grouped$total <- rowSums(combined_dataset_grouped[, c("PPC", "spore spoilage", "no spoilage")])
percentage <- combined_dataset_grouped %>%
  mutate(PPC_percent = PPC / total,
         spore_percent =`spore spoilage` / total) %>%
  select(plantID, PPC_percent, spore_percent)

# Cluster by PPC_percent 
set.seed(1)
kmeans_result_ppc <- kmeans(percentage$PPC_percent, centers = 3)
percentage$PPC_cluster <- kmeans_result_ppc$cluster
percentage$unspoil_percent = 1 - percentage$PPC_percent - percentage$spore_percent
Cluster_1 <- subset(percentage, PPC_cluster == 1)
Cluster_2 <- subset(percentage, PPC_cluster == 2)
Cluster_3 <- subset(percentage, PPC_cluster == 3)
write.csv(percentage,"percentage.csv")
