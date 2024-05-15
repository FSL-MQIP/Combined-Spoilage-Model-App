library(dplyr)
library(tidyr)
library(ggplot2)

test_data = read.csv("test_data_tier_plant.csv")
train_data = read.csv("train_data_tier_plant.csv")
combined_dataset = rbind(test_data,train_data)
combined_dataset = rbind(test_data[,c(5,6,17)],train_data[,c(5,6,17)])

# Not group by Tier 1, 2, 3
# Group by plantID
combined_dataset_grouped <- combined_dataset %>%
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
Cluster_1 <- subset(percentage, PPC_cluster == 1)
Cluster_2 <- subset(percentage, PPC_cluster == 2)
Cluster_3 <- subset(percentage, PPC_cluster == 3)
write.csv(percentage,"percentage.csv")

# Group by Tier 1, 2, 3
Tier_1 = subset(combined_dataset, tier == "1")
Tier_2 = subset(combined_dataset, tier == "2")
Tier_3 = subset(combined_dataset, tier == "3")

Tier_1_pc <- Tier_1 %>%
  group_by(plantID) %>%
  summarize(total_obs = n(),
            ppc_obs = sum(spoilagetype == "PPC"),
            percentage_ppc_spoiled = ppc_obs / total_obs * 100)
Tier_1_pc <- Tier_1_pc %>%
  mutate(tier = "Tier_1")

Tier_2_pc <- Tier_2 %>%
  group_by(plantID) %>%
  summarize(total_obs = n(),
            ppc_obs = sum(spoilagetype == "PPC"),
            percentage_ppc_spoiled = ppc_obs / total_obs * 100)
Tier_2_pc <- Tier_2_pc %>%
  mutate(tier = "Tier_2")

Tier_3_pc <- Tier_3 %>%
  group_by(plantID) %>%
  summarize(total_obs = n(),
            ppc_obs = sum(spoilagetype == "PPC"),
            percentage_ppc_spoiled = ppc_obs / total_obs * 100)
Tier_3_pc <- Tier_3_pc %>%
  mutate(tier = "Tier_3")

result_by_tier = rbind(Tier_1_pc, Tier_2_pc, Tier_3_pc)
write.csv(result_by_tier,"result_by_tier.csv")
