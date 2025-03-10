# Chi-square test
PPC_D7 <- matrix(c(3537, 6463, 44, 223), 
               nrow = 2, 
               byrow = TRUE)
rownames(PPC_D7) <- c("Simulated", "Observed")
colnames(PPC_D7) <- c(">6log", "<6log")
result_PPC_D7 <- chisq.test(PPC_D7)
# p-value < 0.05

PPC_D14 <- matrix(c(8901, 1099, 233, 34), 
                 nrow = 2, 
                 byrow = TRUE)
rownames(PPC_D14) <- c("Simulated", "Observed")
colnames(PPC_D14) <- c(">6log", "<6log")
result_PPC_D14 <- chisq.test(PPC_D14)
# p-value > 0.05

PPC_D21 <- matrix(c(9874, 126, 114, 13), 
                  nrow = 2, 
                  byrow = TRUE)
rownames(PPC_D21) <- c("Simulated", "Observed")
colnames(PPC_D21) <- c(">6log", "<6log")
result_PPC_D21 <- chisq.test(PPC_D21)
# p-value < 0.05

SPC_D7 <- matrix(c(2, 9998, 0, 188), 
                 nrow = 2, 
                 byrow = TRUE)
rownames(SPC_D7) <- c("Simulated", "Observed")
colnames(SPC_D7) <- c(">6log", "<6log")
result_SPC_D7 <- fisher.test(SPC_D7)
# p-value > 0.05

SPC_D14 <- matrix(c(384, 9616, 10, 179), 
                 nrow = 2, 
                 byrow = TRUE)
rownames(SPC_D14) <- c("Simulated", "Observed")
colnames(SPC_D14) <- c(">6log", "<6log")
result_SPC_D14 <- chisq.test(SPC_D14)
# p-value > 0.05

SPC_D21 <- matrix(c(1365, 8635, 39, 94), 
                  nrow = 2, 
                  byrow = TRUE)
rownames(SPC_D21) <- c("Simulated", "Observed")
colnames(SPC_D21) <- c(">6log", "<6log")
result_SPC_D21 <- chisq.test(SPC_D21)
# p-value < 0.05
