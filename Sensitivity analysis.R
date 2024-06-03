library(epiR)

# Sensitivity analysis
Good_Plant_D14_SA = read.csv("Good_Plant_D14_SA.CSV")
SA_result_Good_Plant_D14 = epi.prcc(Good_Plant_D14_SA[,-1], sided.test = 2, conf.level = 0.95)

Good_Plant_D21_SA = read.csv("Good_Plant_D21_SA.CSV")
SA_result_Good_Plant_D21 = epi.prcc(Good_Plant_D21_SA[,-1], sided.test = 2, conf.level = 0.95)

Medium_Plant_D14_SA = read.csv("Medium_Plant_D14_SA.CSV")
SA_result_Medium_Plant_D14 = epi.prcc(Medium_Plant_D14_SA[,-1], sided.test = 2, conf.level = 0.95)

Medium_Plant_D21_SA = read.csv("Medium_Plant_D21_SA.CSV")
SA_result_Medium_Plant_D21 = epi.prcc(Medium_Plant_D21_SA[,-1], sided.test = 2, conf.level = 0.95)

Bad_Plant_D14_SA = read.csv("Bad_Plant_D14_SA.CSV")
SA_result_Bad_Plant_D14 = epi.prcc(Bad_Plant_D14_SA[,-1], sided.test = 2, conf.level = 0.95)

Bad_Plant_D21_SA = read.csv("Bad_Plant_D21_SA.CSV")
SA_result_Bad_Plant_D21 = epi.prcc(Bad_Plant_D21_SA[,-1], sided.test = 2, conf.level = 0.95)

Combined_D14_SA = rbind(Good_Plant_D14_SA, Medium_Plant_D14_SA, Bad_Plant_D14_SA)
SA_result_Combined_D14 = epi.prcc(Combined_D14_SA[,-1], sided.test = 2, conf.level = 0.95)

Combined_D21_SA = rbind(Good_Plant_D21_SA, Medium_Plant_D21_SA, Bad_Plant_D21_SA)
SA_result_Combined_D21 = epi.prcc(Combined_D21_SA[,-1], sided.test = 2, conf.level = 0.95)
