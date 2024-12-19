##################################################
###       Results Empirical Application       ####
##################################################
library(ggplot2)
library(dplyr)
library(stringr)
library(modelconf)
library(GAS)
library(esback)
library(rugarch)
library(xtable)
library(esreg)
library(Rcpp)
library(quantreg)
library(tidyr)
library(ggpubr)
library(murphydiagram)
sourceCpp("scoring_functions.cpp")
source("Function_VaR_VQR.R")
source("GiacominiRossiTest.R")

### Import results
Dates <- lubridate::ymd(read.csv("./Crypto/ETHUSDT_1d.csv", head = TRUE)[-c(1:1001),"OpenTime"])
r_oos <- read.csv("./Crypto/r_oos_ETH_1000.csv", head = FALSE)[,1]
ES_1 <- read.csv("./Crypto/ES1_ETH_1000.csv", head = FALSE)
ES_2 <- read.csv("./Crypto/ES2_ETH_1000.csv", head = FALSE)
VaR_1 <- read.csv("./Crypto/VaR1_ETH_1000.csv", head = FALSE)
VaR_2 <- read.csv("./Crypto/VaR2_ETH_1000.csv", head = FALSE)



# Setup
K <- ncol(ES_1)
a1 <- 0.010
a2 <- 0.025
BackVaRES1 = BackVaRES2 = matrix(0,ncol = 12,nrow = K) 
colnames(BackVaRES1) = colnames(BackVaRES2) = c("Hits", "UC", "CC", "DQ", "VQ", "MFE", "NZ", "ESR_3", "AQL", "AFZG", "ANZ", "AAL")


for (i in 1:K) { 
  print(i)
  set.seed(1234)
  VaRBack1 <- BacktestVaR(r_oos, VaR_1[,i], alpha = a1, Lags = 4)
  VaRBack2 <- BacktestVaR(r_oos, VaR_2[,i], alpha = a2, Lags = 4)

  EBack1 = ESTest(alpha = a1, r_oos, ES_1[,i], VaR_1[,i], conf.level = 0.95,  boot = TRUE, n.boot = 5000)
  EBack2 = ESTest(alpha = a2, r_oos, ES_2[,i], VaR_2[,i], conf.level = 0.95,  boot = TRUE, n.boot = 5000)

  BackVaRES1[i,] = c(mean(r_oos < VaR_1[,i])*100, 
                     VaRBack1$LRuc[2], VaRBack1$LRcc[2],VaRBack1$DQ$pvalue, VaR_VQR(r_oos, VaR_1[,i], a1),
                     EBack1$boot.p.value,
                     cc_backtest(r_oos,  VaR_1[,i], ES_1[,i], alpha  = a1)$pvalue_twosided_simple, 
                     esr_backtest(r_oos, VaR_1[,i], ES_1[,i],alpha  = a1, B = 0, version = 3)$pvalue_onesided_asymptotic,
                     mean(QL(matrix(VaR_1[,i], ncol = 1) ,r_oos, alpha = a1)),
                     mean(FZG(matrix(VaR_1[,i], ncol = 1), matrix(ES_1[,i], ncol = 1), r_oos, alpha = a1)),
                     mean(NZ(matrix(VaR_1[,i], ncol = 1), matrix(ES_1[,i], ncol = 1), r_oos, alpha = a1)),
                     mean(AL(matrix(VaR_1[,i], ncol = 1), matrix(ES_1[,i], ncol = 1), r_oos, alpha = a1)))
  
  BackVaRES2[i,] = c(mean(r_oos < VaR_2[,i])*100, 
                     VaRBack2$LRuc[2], VaRBack2$LRcc[2],VaRBack2$DQ$pvalue, VaR_VQR(r_oos, VaR_2[,i], a2),
                     EBack2$boot.p.value,
                     cc_backtest(r_oos,  VaR_2[,i], ES_2[,i], alpha  = a2)$pvalue_twosided_simple, 
                     esr_backtest(r_oos, VaR_2[,i], ES_2[,i],alpha  = a2, B = 0, version = 3)$pvalue_onesided_asymptotic,
                     mean(QL(matrix(VaR_2[,i], ncol = 1) ,r_oos, alpha = a2)),
                     mean(FZG(matrix(VaR_2[,i], ncol = 1), matrix(ES_2[,i], ncol = 1), r_oos, alpha = a2)),
                     mean(NZ(matrix(VaR_2[,i], ncol = 1), matrix(ES_2[,i], ncol = 1), r_oos, alpha = a2)),
                     mean(AL(matrix(VaR_2[,i], ncol = 1), matrix(ES_2[,i], ncol = 1), r_oos, alpha = a2)))

}

xtable(BackVaRES1, digits = 4)
xtable(BackVaRES2, digits = 4)

# Figure VaR plot

r <- data.frame(t = Dates, tipo = "Returns", values = r_oos)
colnames(VaR_1) <- c("Gray (Normal)", "Gray (Student-t)", "Klaassen (Normal)", "Klaassen (Student-t)", "Haas (Normal)", "Haas (Student-t)")
var <- VaR_1 |> mutate(t = Dates) |> pivot_longer(cols = c("Gray (Normal)", "Gray (Student-t)", "Klaassen (Normal)", "Klaassen (Student-t)", "Haas (Normal)", "Haas (Student-t)"), values_to = "values", names_to = "VaR")
var$VaR <- factor(var$VaR, c("Gray (Normal)", "Gray (Student-t)", "Klaassen (Normal)", "Klaassen (Student-t)", "Haas (Normal)", "Haas (Student-t)"))

ggplot(var) + geom_line(aes(x = t, y = values, color = VaR), linetype = "dashed") + 
  geom_line(data = r, aes(x = t, y = values)) + ylab("Returns") + xlab(" ") + facet_wrap(.~VaR, ncol = 2) + 
  theme_bw() + 
  theme(legend.position = "bottom")

# MCS

pMCS <- 0.25
MCS_type = "t.range"
block_length =  21

### QL
{
MCS_MQL1 = rep(0,ncol(VaR_1))
MQL1 = QL(as.matrix(VaR_1), r_oos, alpha = a1)
colnames(MQL1) = colnames(VaR_1)
aux_MQL1 = estMCS.quick(MQL1, test = MCS_type, B = 5000, l = block_length, alpha = pMCS)
MCS_MQL1[aux_MQL1] = 1

MCS_MQL2 = rep(0,ncol(VaR_2))
MQL2 = QL(as.matrix(VaR_2), r_oos, alpha = a2)
colnames(MQL2) = colnames(VaR_2)
aux_MQL2 = estMCS.quick(MQL2, test = MCS_type, B = 5000, l = block_length, alpha = pMCS)
MCS_MQL2[aux_MQL2] = 1

}

### FZG
{
MCS_MFZG1 = rep(0,ncol(VaR_1))
MFZG1 = FZG(as.matrix(VaR_1), as.matrix(ES_1), r_oos, alpha = a1)
colnames(MFZG1) = colnames(VaR_1)
aux_MFZG1 = estMCS.quick(MFZG1, test = MCS_type, B = 5000, l = block_length, alpha = pMCS)
MCS_MFZG1[aux_MFZG1] = 1

MCS_MFZG2 = rep(0,ncol(VaR_2))
MFZG2 = FZG(as.matrix(VaR_2), as.matrix(ES_2), r_oos, alpha = a2)
colnames(MFZG2) = colnames(VaR_2)
aux_MFZG2 = estMCS.quick(MFZG2, test = MCS_type, B = 5000, l = block_length, alpha = pMCS)
MCS_MFZG2[aux_MFZG2] = 1

}

### NZ
{
  MCS_MNZ1 = rep(0,ncol(VaR_1))
  MNZ1 = FZG(as.matrix(VaR_1), as.matrix(ES_1), r_oos, alpha = a1)
  colnames(MNZ1) = colnames(VaR_1)
  aux_MNZ1 = estMCS.quick(MNZ1, test = MCS_type, B = 5000, l = block_length, alpha = pMCS)
  MCS_MNZ1[aux_MNZ1] = 1
  
  MCS_MNZ2 = rep(0,ncol(VaR_2))
  MNZ2 = FZG(as.matrix(VaR_2), as.matrix(ES_2), r_oos, alpha = a2)
  colnames(MNZ2) = colnames(VaR_2)
  aux_MNZ2 = estMCS.quick(MNZ2, test = MCS_type, B = 5000, l = block_length, alpha = pMCS)
  MCS_MNZ2[aux_MNZ2] = 1
}

### AL
{
  MCS_MAL1 = rep(0,ncol(VaR_1))
  MAL1 = FZG(as.matrix(VaR_1), as.matrix(ES_1), r_oos, alpha = a1)
  colnames(MAL1) = colnames(VaR_1)
  aux_MAL1 = estMCS.quick(MAL1, test = MCS_type, B = 5000, l = block_length, alpha = pMCS)
  MCS_MAL1[aux_MAL1] = 1
  
  MCS_MAL2 = rep(0,ncol(VaR_2))
  MAL2 = FZG(as.matrix(VaR_2), as.matrix(ES_2), r_oos, alpha = a2)
  colnames(MAL2) = colnames(VaR_2)
  aux_MAL2 = estMCS.quick(MAL2, test = MCS_type, B = 5000, l = block_length, alpha = pMCS)
  MCS_MAL2[aux_MAL2] = 1
}


data.frame(QL = c(MCS_MQL1, MCS_MQL2),
           FZG = c(MCS_MFZG1, MCS_MFZG2),
           NZ = c(MCS_MNZ1, MCS_MNZ2),
           AL = c(MCS_MAL1, MCS_MAL2))



# Giacomini and Rossi test
graficos_fluctuations = function(VaR, ES, Ret, Dates, risklevel, fluc_alpha = 0.05, mu_, b, competitors) {
  # b: column of the best model
  m <- round(mu_ * length(Ret))
  days <- Dates[-c(1:(1 + m - 2))]
  a <- ifelse(risklevel == 2.5, 0.025, ifelse(risklevel == 5, 0.050, 0.010))
  p <- ncol(VaR)
  x <- 1:p
  x <- x[competitors]
  LimSup <- c(0, 0, 0)
  LimInf <- c(0, 0, 0)
  Fluctu <- c(0, 0, 0)
  ScoFun <- c("0", "0", "0")
  A_vs_B <- c("0", "0", "0")
  n <- length(days)
  names <- c("Gray (Normal)", "Gray (Student-t)", "Klaassen (Normal)", "Klaassen (Student-t)", "Haas (Normal)", "Haas (Student-t)")

  for (i in x) {
    GR_QL <- fluct_test(QL(matrix(VaR[, b], ncol = 1), Ret, alpha = a), 
                        QL(matrix(VaR[, i], ncol = 1), Ret, alpha = a),
                        mu = mu_, alpha = fluc_alpha, dmv_fullsample = TRUE)
    
    GR_FZG <- fluct_test(FZG(matrix(VaR[, b], ncol = 1), matrix(ES[, b], ncol = 1), Ret, alpha = a), 
                         FZG(matrix(VaR[, i], ncol = 1), matrix(ES[, i], ncol = 1), Ret, alpha = a),
                         mu = mu_, alpha = fluc_alpha, dmv_fullsample = TRUE)
    
    GR_NZ <- fluct_test(NZ(matrix(VaR[, b], ncol = 1), matrix(ES[, b], ncol = 1), Ret, alpha = a), 
                        NZ(matrix(VaR[, i], ncol = 1), matrix(ES[, i], ncol = 1), Ret, alpha = a),
                         mu = mu_, alpha = fluc_alpha, dmv_fullsample = TRUE)
    
    GR_AL <- fluct_test(AL(matrix(VaR[, b], ncol = 1), matrix(ES[, b], ncol = 1), Ret, alpha = a), 
                         AL(matrix(VaR[, i], ncol = 1), matrix(ES[, i], ncol = 1), Ret, alpha = a),
                         mu = mu_, alpha = fluc_alpha, dmv_fullsample = TRUE)
    
  
    
    LimSup <- c(LimSup, rep(GR_QL$cv_sup, n), rep(GR_FZG$cv_sup, n), rep(GR_NZ$cv_sup, n), rep(GR_AL$cv_sup, n))
    LimInf <- c(LimInf, rep(GR_QL$cv_inf, n), rep(GR_FZG$cv_inf, n), rep(GR_NZ$cv_inf, n), rep(GR_AL$cv_inf, n))
    Fluctu <- c(Fluctu, GR_QL$fluc$y, GR_FZG$fluc$y, GR_NZ$fluc$y, GR_AL$fluc$y)
    ScoFun <- c(ScoFun, rep("QL", n), rep("FZG", n), rep("NZ", n), rep("AL", n))
    A_vs_B <- c(A_vs_B, rep(paste0(names[b], " vs. ", names[i]), 4*n))
    
  }
  data_figure <- data.frame(LimSup, LimInf, Fluctu, ScoFun, A_vs_B)
  data_figure <- data_figure[-c(1, 2, 3), ]
  data_figure$days <- rep(days, 4*length(x))
  data_figure$A_vs_B <- factor(data_figure$A_vs_B, unique(data_figure$A_vs_B))
  
  figure <- ggplot(data = data_figure) + 
    #geom_vline(xintercept = as.Date(c("2020-03-12")), color = "red", linetype = "dashed") + 
    geom_line(aes(x = days, Fluctu), color = "green4") + 
    geom_line(aes(x = days, LimSup), color = "black", linetype = "dashed") + 
    geom_line(aes(x = days, LimInf), color = "black", linetype = "dashed") + 
    ylab("Relative performance") + xlab(" ") + facet_grid(ScoFun ~ A_vs_B) +
    theme_bw()
  
  ggsave(filename = paste0("GR_", risklevel), plot = figure, device = "pdf",
         width = 35, height = 21, units = "cm")
  
}

mu_ <- 0.2
b <- 6
competitors <- c(1, 2, 3, 4, 5)
Ret <- r_oos
fluc_alpha <- 0.05

risklevel <- 1
VaR <- VaR_1
ES <- ES_1
graficos_fluctuations(VaR, ES, Ret, Dates, risklevel, fluc_alpha, mu_, b, competitors)

risklevel <- 2
VaR <- VaR_2
ES <- ES_2
graficos_fluctuations(VaR, ES, Ret, Dates, risklevel, fluc_alpha, mu_, b, competitors)



