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
library(moments)
library(xtable)
library(gridExtra)
library(grid)
library(scales)

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


###################################
###    Table 1 and Figure 1     ###    
###################################

BTC = read.csv("./Data/BTCUSDT_1d.csv") %>% 
  mutate(date = as.Date(OpenTime)) %>% arrange(date) %>% mutate(ret = c(NA,diff(log(Close))*100)) %>% 
  dplyr::select(date, ret) %>% rename(BTC = ret) |> drop_na()

ETH = read.csv("./Data/ETHUSDT_1d.csv") %>% 
  mutate(date = as.Date(OpenTime)) %>% arrange(date) %>% mutate(ret = c(NA,diff(log(Close))*100)) %>% 
  dplyr::select(date, ret) %>% rename(ETH = ret) |> drop_na()

Crypto = BTC %>% left_join(ETH, by = "date") %>% dplyr::select(date, BTC, ETH) 

Descriptive = function(ret){
  c(min(ret, na.rm = TRUE), quantile(ret, 0.25, na.rm = TRUE), quantile(ret, 0.5, na.rm = TRUE), 
    mean(ret, na.rm = TRUE), quantile(ret, 0.75, na.rm = TRUE), max(ret, na.rm = TRUE),
    sd(ret, na.rm = TRUE), skewness(ret, na.rm = TRUE), kurtosis(ret, na.rm = TRUE))
}
N = dim(Crypto)[2] - 1 
Results = matrix(0, ncol = 9, nrow = N)
for (i in 1:N) {
  Results[i, ] = Descriptive(Crypto[, i + 1])
}

colnames(Results) = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "Sd", "Skew", "Kurtosis")
row.names(Results) = colnames(Crypto[, -1])
tabela = xtable(Results,caption = "Descriptive statistics of daily returns", digits = 2)

InS <- 1000

{
  
  p1_1 = ggplot(BTC, aes(date, BTC)) + geom_line(colour = "green4") +
    xlab(" ") + ylab("Bitcoin") + geom_vline(xintercept = Crypto[InS + 1,"date"], linetype = "dashed") + 
    scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
    theme_bw() + theme(axis.text.x = element_text(size = 9,face = "bold")) 
  
  p1_2 = ggplot(ETH, aes(date, ETH)) + geom_line(colour = "green4") +
    xlab(" ") + ylab("Ethereum") + geom_vline(xintercept = Crypto[InS + 1,"date"], linetype = "dashed") + 
    scale_x_date(date_breaks = "1 year", date_labels =  "%Y") +
    theme_bw() + theme(axis.text.x = element_text(size = 9,face = "bold"))
}

######### Autocorrelations Generalised Barlett
{
  
  gamma <- function(x,h){
    n <- length(x)
    h <- abs(h)
    x <- x - mean(x)
    gamma <- sum(x[1:(n - h)]*x[(h + 1):n])/n
  }
  rho <- function(x,h) rho <- gamma(x,h)/gamma(x,0)
  
  
  x = BTC$BTC
  n = length(x)
  nlag <- 50 
  acf.val <- sapply(c(1:nlag),function(h) rho(x,h))
  x2 <- x^2
  var <- 1 + (sapply(c(1:nlag),function(h) gamma(x2,h)))/gamma(x,0)^2
  band <- sqrt(var/n)
  
  conf.level <- 0.95
  ciline <- qnorm((1 - conf.level)/2)/sqrt(length(x))
  bacf <- acf(x, plot = FALSE, 50)  
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  p2_1 <- ggplot(data = bacfdf[-1,], mapping = aes(x = lag, y = acf)) +
    geom_bar(stat = "identity", position = "identity", fill = "green4") + ylab("Bitcoin") + 
    ylim(c(-0.15,0.15)) +
    geom_line(aes(y = -1.96*band, x = 1:50)) +
    geom_line(aes(y = -1.96/sqrt(n), x = 1:50), linetype = "dotted") +
    geom_line(aes(y = 1.96/sqrt(n), x = 1:50), linetype = "dotted") +
    geom_line(aes(y = 1.96*band, x = 1:50)) + 
    theme_bw() + 
    theme(legend.position = "none")
  
  x = ETH$ETH
  n = length(x)
  nlag <- 50 
  acf.val <- sapply(c(1:nlag),function(h) rho(x,h))
  x2 <- x^2
  var <- 1 + (sapply(c(1:nlag),function(h) gamma(x2,h)))/gamma(x,0)^2
  band <- sqrt(var/n)
  
  conf.level <- 0.95
  ciline <- qnorm((1 - conf.level)/2)/sqrt(length(x))
  bacf <- acf(x, plot = FALSE, 50)  
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  p2_2 <- ggplot(data = bacfdf[-1,], mapping = aes(x = lag, y = acf)) +
    geom_bar(stat = "identity", position = "identity", fill = "green4") + ylab("Ethereum") + 
    ylim(c(-0.15,0.15)) +
    geom_line(aes(y = -1.96*band, x = 1:50)) +
    geom_line(aes(y = -1.96/sqrt(n), x = 1:50), linetype = "dotted") +
    geom_line(aes(y = 1.96/sqrt(n), x = 1:50), linetype = "dotted") +
    geom_line(aes(y = 1.96*band, x = 1:50)) + 
    theme_bw() + 
    theme(legend.position = "none")
}

######### Autocorrelations Barlett
{
  retornos = BTC$BTC
  inventories = retornos^2
  n <- length(inventories)
  mean.inventories <- sum(inventories)/n
  # Express the data in deviations from the mean
  z.bar <- rep(mean.inventories,n)
  deviations <- inventories - z.bar
  # Calculate the sum of squared deviations from the mean
  squaredDeviations <- deviations^2
  sumOfSquaredDeviations <- sum(squaredDeviations)
  # Create empty vector to store autocorrelation coefficients
  r <- c()
  # Use a for loop to fill the vector with the coefficients
  for (k in 1:n) {
    ends <- n - k
    starts <- 1 + k
    r[k] <- sum(deviations[1:(ends)]*deviations[(starts):(n)])/sumOfSquaredDeviations
  }
  # Create empty vector to store Bartlett's standard errors
  bart.error <- c()
  # Use a for loop to fill the vector with the standard errors
  for (k in 1:n) {
    ends <- k - 1
    bart.error[k] <- ((1 + sum((2*r[0:(ends)]^2)))^0.5)*(n^-0.5)
  }
  
  bacf <- acf(retornos^2, plot = FALSE, 50)  
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  
  
  p3_1 <- ggplot(data = bacfdf[-1,], mapping = aes(x = lag, y = acf)) +
    geom_bar(stat = "identity", position = "identity", fill = "green4") + ylab(expression(Bitcoin^{2})) + 
    geom_line(aes(y = -1.96*bart.error[1:50], x = 1:50)) +
    geom_line(aes(y = 1.96*bart.error[1:50], x = 1:50)) + theme_bw() + 
    theme(legend.position = "none")
  
  retornos = ETH$ETH
  inventories = retornos^2
  n <- length(inventories)
  mean.inventories <- sum(inventories)/n
  # Express the data in deviations from the mean
  z.bar <- rep(mean.inventories,n)
  deviations <- inventories - z.bar
  # Calculate the sum of squared deviations from the mean
  squaredDeviations <- deviations^2
  sumOfSquaredDeviations <- sum(squaredDeviations)
  # Create empty vector to store autocorrelation coefficients
  r <- c()
  # Use a for loop to fill the vector with the coefficients
  for (k in 1:n) {
    ends <- n - k
    starts <- 1 + k
    r[k] <- sum(deviations[1:(ends)]*deviations[(starts):(n)])/sumOfSquaredDeviations
  }
  # Create empty vector to store Bartlett's standard errors
  bart.error <- c()
  # Use a for loop to fill the vector with the standard errors
  for (k in 1:n) {
    ends <- k - 1
    bart.error[k] <- ((1 + sum((2*r[0:(ends)]^2)))^0.5)*(n^-0.5)
  }
  
  bacf <- acf(retornos^2, plot = FALSE, 50)  
  bacfdf <- with(bacf, data.frame(lag, acf))
  
  
  
  p3_2 <- ggplot(data = bacfdf[-1,], mapping = aes(x = lag, y = acf)) +
    geom_bar(stat = "identity", position = "identity", fill = "green4") + ylab(expression(Ethereum^{2})) + 
    geom_line(aes(y = -1.96*bart.error[1:50], x = 1:50)) +
    geom_line(aes(y = 1.96*bart.error[1:50], x = 1:50)) + theme_bw() + 
    theme(legend.position = "none")
  
}

pdf("crypto_figures.pdf", paper = "a4r", width = 18, height = 8) 
pushViewport(viewport(layout = grid.layout(2, 3)))
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
print(p1_1, vp = vplayout(1, 1))
print(p2_1, vp = vplayout(1, 2))
print(p3_1, vp = vplayout(1, 3))
print(p1_2, vp = vplayout(2, 1))
print(p2_2, vp = vplayout(2, 2))
print(p3_2, vp = vplayout(2, 3))
dev.off()
