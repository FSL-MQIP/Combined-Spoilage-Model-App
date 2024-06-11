library(epiR)
library(gridExtra)

# Sensitivity analysis
Good_Plant_D14_SA = read.csv("Good_Plant_D14_SA.CSV")
SA_result_Good_Plant_D14 = epi.prcc(Good_Plant_D14_SA[,-1], sided.test = 2, conf.level = 0.95)
SA_result_Good_Plant_D14$var = c("t_F","T_F","t_T","T_T","t_S","T_S","t_T2","T_T2","T_H","P_ppc","N0")
SA_result_Good_Plant_D14$abs_PRCC_value = abs(SA_result_Good_Plant_D14$est)

Good_Plant_D21_SA = read.csv("Good_Plant_D21_SA.CSV")
SA_result_Good_Plant_D21 = epi.prcc(Good_Plant_D21_SA[,-1], sided.test = 2, conf.level = 0.95)
SA_result_Good_Plant_D21$var = c("t_F","T_F","t_T","T_T","t_S","T_S","t_T2","T_T2","T_H","P_ppc","N0")
SA_result_Good_Plant_D21$abs_PRCC_value = abs(SA_result_Good_Plant_D21$est)

Medium_Plant_D14_SA = read.csv("Medium_Plant_D14_SA.CSV")
SA_result_Medium_Plant_D14 = epi.prcc(Medium_Plant_D14_SA[,-1], sided.test = 2, conf.level = 0.95)
SA_result_Medium_Plant_D14$var = c("t_F","T_F","t_T","T_T","t_S","T_S","t_T2","T_T2","T_H","P_ppc","N0")
SA_result_Medium_Plant_D14$abs_PRCC_value = abs(SA_result_Medium_Plant_D14$est)

Medium_Plant_D21_SA = read.csv("Medium_Plant_D21_SA.CSV")
SA_result_Medium_Plant_D21 = epi.prcc(Medium_Plant_D21_SA[,-1], sided.test = 2, conf.level = 0.95)
SA_result_Medium_Plant_D21$var = c("t_F","T_F","t_T","T_T","t_S","T_S","t_T2","T_T2","T_H","P_ppc","N0")
SA_result_Medium_Plant_D21$abs_PRCC_value = abs(SA_result_Medium_Plant_D21$est)

Bad_Plant_D14_SA = read.csv("Bad_Plant_D14_SA.CSV")
SA_result_Bad_Plant_D14 = epi.prcc(Bad_Plant_D14_SA[,-1], sided.test = 2, conf.level = 0.95)
SA_result_Bad_Plant_D14$var = c("t_F","T_F","t_T","T_T","t_S","T_S","t_T2","T_T2","T_H","P_ppc","N0")
SA_result_Bad_Plant_D14$abs_PRCC_value = abs(SA_result_Bad_Plant_D14$est)

Bad_Plant_D21_SA = read.csv("Bad_Plant_D21_SA.CSV")
SA_result_Bad_Plant_D21 = epi.prcc(Bad_Plant_D21_SA[,-1], sided.test = 2, conf.level = 0.95)
SA_result_Bad_Plant_D21$var = c("t_F","T_F","t_T","T_T","t_S","T_S","t_T2","T_T2","T_H","P_ppc","N0")
SA_result_Bad_Plant_D21$abs_PRCC_value = abs(SA_result_Bad_Plant_D21$est)

# Tornado plot 
SA_result_Good_Plant_D21_filtered = SA_result_Good_Plant_D21[SA_result_Good_Plant_D21$p.value < 0.05, ]
SA_result_Medium_Plant_D21_filtered = SA_result_Medium_Plant_D21[SA_result_Medium_Plant_D21$p.value < 0.05, ]
SA_result_Bad_Plant_D21_filtered = SA_result_Bad_Plant_D21[SA_result_Bad_Plant_D21$p.value < 0.05, ]

SA_result_Bad_Plant_D21_filtered$abs_prcc = abs(SA_result_Bad_Plant_D21_filtered$est)
SA_result_Bad_Plant_D21_filtered = SA_result_Bad_Plant_D21_filtered[order(-SA_result_Bad_Plant_D21_filtered$abs_prcc),]
SA_result_Bad_Plant_D21_filtered$var = c("Initial contamination level", 
                                         "Consumer home storage temperature", 
                                         "Transportation temperature between processing facility and retail", 
                                         "Percent spoiled by PPC", 
                                         "Retail storage temperature",
                                         "Retail storage time", 
                                         "Transportation time between processing facility and retail"
                                          )


plot_Bad = ggplot(SA_result_Bad_Plant_D21_filtered, aes(x = est, y = reorder(var, abs_prcc))) +
  geom_bar(stat = "identity", fill = "grey") +
  labs(x = "Partial Rank Correlation Coefficient", y = "Parameters", title = "Bad Plant") +
  theme_minimal() + 
  xlim(-0.1,0.6)

grid.arrange(plot_Good, plot_Medium, plot_Bad, ncol = 1)
