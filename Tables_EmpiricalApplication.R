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
sourceCpp("Aux/scoring_functions.cpp")
source("Aux/Function_VaR_VQR.R")

### Import results
### Use ETH or BTC
Dates <- lubridate::ymd(read.csv("./Data/ETHUSDT_1d.csv", head = TRUE)[-c(1:1001),"OpenTime"])
r_oos <- read.csv("r_oos_ETH_1000.csv", head = FALSE)[,1]
ES_1 <- read.csv("ES1_ETH_1000.csv", head = FALSE)
ES_2 <- read.csv("ES2_ETH_1000.csv", head = FALSE)
VaR_1 <- read.csv("VaR1_ETH_1000.csv", head = FALSE)
VaR_2 <- read.csv("VaR2_ETH_1000.csv", head = FALSE)



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

# 1: Belongs to the MCS
# 0: Does not belong to the MCS