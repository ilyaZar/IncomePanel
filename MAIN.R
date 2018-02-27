#### Main Kneib Implementation ####
setwd("/home/chief/Dropbox/Projects/Research/Basti1/Kneib/CrossSection")
# rm(list = ls())
library(MASS)
library(spdep)
library(mvtnorm)
library(Matrix)
library(GB2)
library(optimx)
library(matrixcalc)
library(plyr)
library(coda)
library(ggplot2)
library(ggmcmc)
library(gridExtra)
source("HelperSamplers.R")
#
#
#
#
#
set.seed(123)
n      <- 100000
kT     <- 10
Ncross <- 10
time   <- 3

tauAT <- 1
tauBT <- 1
tauPT <- 1
tauQT <- 1

dimBA <- 2
dimBB <- 2
dimBP <- 2
dimBQ <- 2
dimB  <- dimBA + dimBB + dimBP + dimBQ

dimCA <- 10
dimCB <- 1
dimCP <- 1
dimCQ <- 1

betaAT <- as.vector(rmvnorm(1, mean = rep(0, times = dimBA),
                            sigma = diag(x = tauAT^2,
                            nrow = dimBA, ncol = dimBA)))
betaBT <- as.vector(rmvnorm(1, mean = rep(0, times = dimBB),
                            sigma = diag(x = tauBT^2,
                            nrow = dimBB, ncol = dimBB)))
betaPT <- as.vector(rmvnorm(1, mean = rep(0, times = dimBP),
                            sigma = diag(x = tauPT^2,
                            nrow = dimBP, ncol = dimBP)))
betaQT <- as.vector(rmvnorm(1, mean = rep(0, times = dimBQ),
                            sigma = diag(x = tauQT^2,
                            nrow = dimBQ, ncol = dimBQ)))
#(1/rgamma(n = 1, shape = 0.01, scale = 1/(0.01)))^(-0.5)
# aT    <- 1/rgamma(n = 1, shape = tauAT, scale = 1/phiAT)
# ceptAT <- (1:dimCA)*0.1
ceptBT <- -((1:dimCB)*0.1)
ceptPT <- -((1:dimCP)*0.3)
ceptQT <- (1:dimCQ)*0.3
ceptAT <- rnorm(dimCA, 0, 5)
# ceptBT <- rep(-0.5, times = dimCB)
# ceptPT <- rep(-0.5, times = dimCP)
# ceptQT <- rep(0.5,  times = dimCQ)
# indiAT <- 1:Ncross
# indiBT <- -(1:Ncross/2)
# indiPT <- 1:Ncross/3
# indiQT <- -(1:Ncross/4)

trueVAL <- c(betaAT, betaBT, betaPT, betaQT, ceptAT, ceptBT, ceptPT, ceptQT)

aT    <- 1.5
bT    <- 150
pT    <- 2.5
qT    <- 3.5
SDreg <- 0.1 # Up to 0.7 with one regressor 

if (dimBA > 1) {
  XA  <- matrix(0, nrow = Ncross*time, ncol = dimBA)
  XA[, 1:(dimBA - 1)] <- matrix(rnorm(n = Ncross*time*(dimBA - 1), 0, 1), nrow = Ncross*time, ncol = (dimBA - 1))
  regMean <- log(aT) - ceptAT 
  regMean <- (regMean - as.matrix(XA[, 1:(dimBA - 1)]) %*% betaAT[1:(dimBA - 1)])
  regMean <- regMean/betaAT[dimBA]
  XA[, dimBA] <- rnorm(n = Ncross*time, mean = regMean, sd = SDreg)
} else {
  XA  <- matrix(rnorm(n = Ncross*time*dimBA, mean = (log(aT) - ceptAT)/betaAT, sd = SDreg),
                nrow = Ncross*time, ncol = dimBA)
}

if (dimBB > 1) {
  XB  <- matrix(0, nrow = Ncross*time, ncol = dimBB)
  XB[, 1:(dimBB - 1)] <- matrix(rnorm(n = Ncross*time*(dimBB - 1), 0, 1), nrow = Ncross*time, ncol = (dimBB - 1))
  regMean <- log(bT) - ceptBT
  regMean <- (regMean - as.matrix(XB[, 1:(dimBB - 1)]) %*% betaBT[1:(dimBB - 1)])
  regMean <- regMean/betaBT[dimBB]
  XB[, dimBB] <- rnorm(n = Ncross*time, mean = regMean, sd = SDreg)
} else {
  XB  <- matrix(rnorm(n = Ncross*time*dimBB, mean = (log(bT) - ceptBT)/betaBT, sd = SDreg),
                nrow = Ncross*time, ncol = dimBB)
}

if (dimBP > 1) {
  XP  <- matrix(0, nrow = Ncross*time, ncol = dimBP)
  XP[, 1:(dimBP - 1)] <- matrix(rnorm(n = Ncross*time*(dimBP - 1), 0, 1), nrow = Ncross*time, ncol = (dimBP - 1))
  regMean <- log(pT) - ceptPT
  regMean <- (regMean - as.matrix(XP[, 1:(dimBP - 1)]) %*% betaPT[1:(dimBP - 1)])
  regMean <- regMean/betaPT[dimBP]
  XP[, dimBP] <- rnorm(n = Ncross*time, mean = regMean, sd = SDreg)
} else {
  XP  <- matrix(rnorm(n = Ncross*time*dimBP, mean = (log(pT) - ceptPT)/betaPT, sd = SDreg),
                nrow = Ncross*time, ncol = dimBP)
}

if (dimBQ > 1) {
  XQ  <- matrix(0, nrow = Ncross*time, ncol = dimBQ)
  XQ[, 1:(dimBQ - 1)] <- matrix(rnorm(n = Ncross*time*(dimBQ - 1), 0, 1), nrow = Ncross*time, ncol = (dimBQ - 1))
  regMean <- log(qT) - ceptQT
  regMean <- (regMean - as.matrix(XQ[, 1:(dimBQ - 1)]) %*% betaQT[1:(dimBQ - 1)])
  regMean <- regMean/betaQT[dimBQ]
  XQ[, dimBQ] <- rnorm(n = Ncross*time, mean = regMean, sd = SDreg)
} else {
  XQ  <- matrix(rnorm(n = Ncross*time*dimBQ, mean = (log(qT) - ceptQT)/betaQT, sd = SDreg),
                nrow = Ncross*time, ncol = dimBQ)
}

aT    <- exp(ceptAT + XA %*% betaAT)
bT    <- exp(ceptBT + XB %*% betaBT)
pT    <- exp(ceptPT + XP %*% betaPT)
qT    <- exp(ceptQT + XQ %*% betaQT)
# 
# 
# 
# 
# Defining print strings:
ARSTRING     <- " Acceptance rate MH:     %4f \n"
PARSTRING1   <- paste(" Parameter    :  %s     \n")
PARSTRING2   <- paste(rep(c("betaA ", "betaB ", "betaP ", "betaQ ", "ceptBA ",
                          "betaBB ", "betaBP ", "ceptBQ "),
                        times = c(dimBA, dimBB, dimBP, dimBQ, dimCA, dimCB, dimCP, dimCQ)), collapse = "")
TRUESTRING1  <- paste(" True values  :   %s     \n")
TRUESTRING2  <- paste(trueVAL, collapse = "      ")
CURSTRING    <- paste(" Current state:",
                     paste(rep("%.4f  ", times = length(trueVAL)), collapse = ""),
                     "\n",
                     collapse = "")
MEANSTRING   <- paste(" Current mean :",
                     paste(rep("%.4f  ", times = length(trueVAL)), collapse = ""),
                     "\n \n",
                     collapse = "")
#
#
#
#
#
# Defining index references for "Samples"-type storing matrices
indBA <- 1:dimBA
indBB <- (dimBA + 1):(dimBA +  dimBB)
indBP <- (dimBA +  dimBB + 1):(dimBA +  dimBB + dimBP)
indBQ <- (dimBA +  dimBB + dimBP + 1):(dimBA +  dimBB + dimBP + dimBQ)
indCA <- (dimBA +  dimBB + dimBP + dimBQ) + 1:dimCA
indCB <- (dimBA +  dimBB + dimBP + dimBQ) + dimCA + 1:dimCB
indCP <- (dimBA +  dimBB + dimBP + dimBQ) + dimCA + dimCB + 1:dimCP
indCQ <- (dimBA +  dimBB + dimBP + dimBQ) + dimCA + dimCB + dimCP + 1:dimCQ
#
#
#
#
# Data Simulation ==============================================================
set.seed(123)
dataSIM <- list()
zT      <- list()
for (j in 1:time) {
  for (i in 1:Ncross) {
    ind <- (j - 1)*Ncross + i
    z <- qgb2(prob = seq(from = 0, to = 1 - (1/kT), length.out = kT),
              shape1 = aT[ind], scale = bT[ind],
              shape2 = pT[ind], shape3 = qT[ind])
    z <- c(z, Inf)
    zT[[ind]] <- z
    dataSIM[[ind]] <- DGP2(N = n, A = aT[ind], B = bT[ind],
                           P = pT[ind], Q = qT[ind], Z = zT[[ind]])
    # Nbar <- dataSIM[[ind]]
    # out <- optimx(par = trueVAL[1:4], fn = lpost, gr = NULL, hess = NULL,
    #               itnmax = 5000, hessian = TRUE, method = "BFGS",
    #               nbar = Nbar, Z = zT, k = kT, 
    #               Xa = XA[ind], Xb = XB[ind], Xp = XP[ind], Xq = XQ[ind],
    #               control = list(maximize = TRUE,
    #                             reltol  = 10^(-12),
    #                             # kkttol  = .Machine$double.eps^(1/2),
    #                             kkt2tol = 10^(-8))
    # )
    # VCMest  <- -solve(attr(out, "details")[1,]$nhatend)
    # print(out)
    # print(VCMest)
  }
}
#
#
#
#
#
# Estimating proposal VCM
# lpost(nbars = dataSIM, pars = trueVAL, Z = zT, k = kT,
#       Xa = XA, Xb = XB, Xp = XP, Xq = XQ, N = Ncross, t = time,
#       indba = indBA, indbb = indBB, indbp = indBP, indbq = indBQ,
#       indca = indCA, indcb = indCB, indcp = indCP, indcq = indCQ)
out <- optimx(par = trueVAL, fn = lpost, gr = NULL, hess = NULL,
              itnmax = 5000, hessian = TRUE, method = "BFGS",
              nbar = dataSIM, Z = zT, k = kT, t = time,
              Xa = XA, Xb = XB, Xp = XP, Xq = XQ, N = Ncross,
              indba = indBA, indbb = indBB, indbp = indBP, indbq = indBQ,
              indca = indCA, indcb = indCB, indcp = indCP, indcq = indCQ,
              control = list(maximize = TRUE,
                             reltol  = 10^(-12),
                             # kkttol  = .Machine$double.eps^(1/2),
                             kkt2tol = 10^(-8)
              )
)
View(out)
print(trueVAL)
VCMest  <- -solve(attr(out, "details")[1,]$nhatend)
Nsims <- 1
SamplesBeta <- list()
SamplesCept <- list()
for (i in 1:Nsims) {
  print(i)
  # MCMC Bookkepping & starting values:
  M         <- 10000
  Burn      <- 1# 25000# round(M/2)
  aR        <- numeric(M)
  ARall     <- numeric(M)
  # Starting Values:
  updownInd  <- sample(rep(c(1, -1), times = length(trueVAL)))
  PercERR <- 0
  betaA <- betaAT + (updownInd[indBA] * abs(betaAT*PercERR))
  betaB <- betaBT + (updownInd[indBB] * abs(betaBT*PercERR))
  betaP <- betaPT + (updownInd[indBP] * abs(betaPT*PercERR))
  betaQ <- betaQT + (updownInd[indBQ] * abs(betaQT*PercERR))
  ceptA <- ceptAT + (updownInd[indCA] * abs(ceptAT*PercERR))
  ceptB <- ceptBT + (updownInd[indCB] * abs(ceptBT*PercERR))
  ceptP <- ceptPT + (updownInd[indCP] * abs(ceptPT*PercERR))
  ceptQ <- ceptQT + (updownInd[indCQ] * abs(ceptQT*PercERR))
  
  perturbVAL <- c(betaA, betaB, betaP, betaQ, ceptA, ceptB, ceptP, ceptQ)
  
  SamplesBeta[[i]]      <- matrix(0, nrow = M, ncol = dimBA + dimBB + dimBP + dimBQ)
  SamplesCept[[i]]      <- matrix(0, nrow = M, ncol = dimCA + dimCB + dimCP + dimCQ)
  SamplesBeta[[i]][1, ] <- c(betaA, betaB, betaP, betaQ) 
  SamplesCept[[i]][1, ] <- c(ceptA, ceptB, ceptP, ceptQ)
  #
  #
  #
  #
  #
  ### MCMC: ####
  for (m in 2:M) {
    # Pre-comute
    # if (m >= 2000) {
    #   VCMest <- cov(cbind(SamplesBeta[1:(m - 1), ], SamplesCept[1:(m - 1), ]))
    # }
    # 
    # 
    # 
    # 
    #
    betaRES <- betaMHall(nbars = dataSIM,
                         ceptAC = ceptA, ceptBC = ceptB, 
                         ceptPC = ceptP, ceptQC = ceptQ,
                         betaAC = betaA, betaBC = betaB,
                         betaPC = betaP, betaQC = betaQ,
                         Xa = XA, Xb = XB, Xp = XP, Xq = XQ, 
                         N = Ncross, Z = zT, k = kT, VCM = VCMest, t = time, 
                         indba = indBA, indbb = indBB, indbp = indBP, indbq = indBQ,
                         indca = indCA, indcb = indCB, indcp = indCP, indcq = indCQ)
    betaA   <- betaRES[[1]][indBA]
    betaB   <- betaRES[[1]][indBB]
    betaP   <- betaRES[[1]][indBP]
    betaQ   <- betaRES[[1]][indBQ]
    ceptA   <- betaRES[[1]][indCA]
    ceptB   <- betaRES[[1]][indCB]
    ceptP   <- betaRES[[1]][indCP]
    ceptQ   <- betaRES[[1]][indCQ]
    aR[m]   <- betaRES[[2]]
    ARall[m] <- mean(aR[1:m])
    #
    #
    #
    #
    #
    SamplesBeta[[i]][m, indBA] <- betaA
    SamplesBeta[[i]][m, indBB] <- betaB
    SamplesBeta[[i]][m, indBP] <- betaP
    SamplesBeta[[i]][m, indBQ] <- betaQ
    SamplesCept[[i]][m, indCA - dimB] <- ceptA
    SamplesCept[[i]][m, indCB - dimB] <- ceptB
    SamplesCept[[i]][m, indCP - dimB] <- ceptP
    SamplesCept[[i]][m, indCQ - dimB] <- ceptQ
    #
    #
    #
    #
    #
    print(m)
    # cat(sprintf(
    #   "################################################################################\n"))
    # cat(sprintf(" Iteration:     %d/%d completed.\n \n",
    #             m, M))
    # cat(sprintf(PARSTRING1, PARSTRING2))
    # cat(sprintf(TRUESTRING1, TRUESTRING2))
    # cat(sprintf(CURSTRING,
    #             SamplesBeta[m, indBA], SamplesBeta[m, indBB],
    #             SamplesBeta[m, indBP], SamplesBeta[m, indBQ],
    #             SamplesCept[m, indCA - dimB], SamplesCept[m, indCB - dimB],
    #             SamplesCept[m, indCP - dimB], SamplesCept[m, indCQ - dimB]
    # )
    # )
    # # cat(sprintf(" Proposed next state of the Markov chain: %.4f %.4f %.4f \n",
    # #             thp[kk, 1], thp[kk, 2], thp[kk, 3]))
    # #             posterior # of the Markov chain
    # #             #
    # cat(sprintf(MEANSTRING,
    #             mean(SamplesBeta[0:m, indBA]), mean(SamplesBeta[0:m, indBB]),
    #             mean(SamplesBeta[0:m, indBP]), mean(SamplesBeta[0:m, indBQ]),
    #             mean(SamplesCept[0:m, indCA - dimB]), mean(SamplesCept[0:m, indCB - dimB]),
    #             mean(SamplesCept[0:m, indCP - dimB]), mean(SamplesCept[0:m, indCQ - dimB])
    # )
    # )
    # cat(sprintf(ARSTRING,
    #             ARall[m]
    #             )
    #     )
    # cat(sprintf(" acceptance rate rho:    %.4f \n",
    #             aRRho[m]))
    # cat(sprintf(" acceptance rate theta:  %.4f \n",
    #             aRTheta[m]))
    # cat(sprintf(
    #   "################################################################################\n"))
  }
}
sd
# View(SamplesBeta[[1]])
resMCMC <- cbind(SamplesBeta[[1]][,1:8])
colnames(resMCMC) <- c("betaA1","betaA2","betaB1","betaB2",
                       "betaP1","betaP2","betaQ1","betaQ2")
trueTable <- as.matrix(trueVAL[1:8])
colnames(trueTable) <- "True Val."
resMCMC <- mcmc(resMCMC)
effectiveSize(resMCMC)
MCMCtable <- xtable(cbind(trueTable,
                          summary(resMCMC)[[1]][,-(3:4)],
                          summary(resMCMC)[[2]][,-(2:4)]),
                    auto = TRUE, digits = 4)
MCMCtable
# resMCMC <- as.mcmc.list(resMCMC)
# testMCMC <- ggs(resMCMC)
# ggs_autocorrelation(testMCMC, family = "betaA")
# ggs_density(testMCMC, family = "betaA")
# f1 <- ggs_traceplot(ggs(resMCMC, family = "betaA")) + 
#   theme_fivethirtyeight()
# f2 <- ggs_density(ggs(resMCMC, family = "betaA")) + 
#   theme_solarized(light = TRUE)
# f3 <- ggs_autocorrelation(ggs(resMCMC, family = "betaA")) + 
#   theme_solarized(light = TRUE)
# grid.arrange(f2, f1, f3, ncol = 3, nrow = 2)
# autocorr.plot(resMCMC)
#  crosscorr.plot(resMCMC)
#  
#  
#  
#  
# 
estimateSims <- list()
for (i in 1:Nsims) {
  estimateSims[[i]] <- apply(SamplesBeta[[i]], 2, mean)
}
NumErrMat <- matrix(unlist(estimateSims), ncol = 8, byrow = TRUE)
apply(NumErrMat, 2, sd)
#
SamplesBetaBigSim <- as.matrix(SamplesBeta[[1]])
SamplesCeptBigSim <- as.matrix(SamplesCept[[1]])
SamplesBeta <- SamplesBeta[[1]]
SamplesCept <- SamplesCept[[1]]
layout(matrix(c(1,1,2,3), 2, 2, byrow = F))
title <- list()
title[[1]] <- expression(paste("histogram ", beta,A))
title[[3]] <- expression(paste("histogram ", beta,B))
title[[5]] <- expression(paste("histogram ", beta,P))
title[[7]] <- expression(paste("histogram ", beta,Q))
title[[2]] <- expression(paste("histogram ", beta,A))
title[[4]] <- expression(paste("histogram ", beta,B))
title[[6]] <- expression(paste("histogram ", beta,P))
title[[8]] <- expression(paste("histogram ", beta,Q))
for (i in 1:8) {
  hist(SamplesBeta[Burn:M, i], main = "", xlab = "mean (blue)      true value (red)")
  abline(v = trueVAL[i], col = "red")
  abline(v = mean(SamplesBeta[Burn:M, i]), col = "blue")
  plot(SamplesBeta[1:M, i], type = "l", ylab = "value",
       main  = "full trace", xlab = "# iteration")
  plot(SamplesBeta[Burn:M, i], type = "l", ylab = "value",
       main = "trace after burn-in", xlab = "# iteration")
}
# title[[1]] <- expression(paste("histogram ", intercept, beta,A))
# title[[2]] <- expression(paste("histogram ", intercept, beta,B))
# title[[3]] <- expression(paste("histogram ", intercept, beta,P))
# title[[4]] <- expression(paste("histogram ", intercept, beta,Q))
for (i in 1:4) {
  hist(SamplesCept[Burn:M, i], main = "", xlab = "mean (blue)      true value (red)")
  abline(v = trueVAL[i + 8], col = "red")
  abline(v = mean(SamplesCept[Burn:M, i]), col = "blue")
  plot(SamplesCept[1:M, i], type = "l", ylab = "value",
       main  = "full trace", xlab = "# iteration")
  plot(SamplesCept[Burn:M, i], type = "l", ylab = "value",
       main = "trace after burn-in", xlab = "# iteration")
}
# confLOW <- function(vec){
#   sortVec <- sort(vec)
#   minVal  <- sortVec[round(0.025*length(sortVec), digits = 0)] #round(0.025*length(sortVec), digits = 0)
#   return(minVal)
# }
# confUPP <- function(vec){
#   sortVec <- sort(vec)
#   maxVal  <- sortVec[round(0.975*length(sortVec), digits = 0)]
#   return(maxVal)
# }
# round(cbind(trueVAL[1:8],
#       as.vector(apply(SamplesBeta[Burn:M,], 2, mean)),
#       as.vector(apply(SamplesBeta[Burn:M,], 2, sd)),
#       as.vector(apply(SamplesBeta[Burn:M,], 2, confLOW)),
#       as.vector(apply(SamplesBeta[Burn:M,], 2, confUPP))
#       ), digits = 4)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
# betaRES <- betaAMH(nbar = Nbar,
#                   betaAC = betaA, betaBC = betaB,
#                   betaPC = betaP, betaQC = betaQ,
#                   Xa = XA, Xb = XB, Xp = XP, Xq = XQ,
#                   Z = zT, k = kT, VCM = VCMestA)
# betaA      <- betaRES[[1]][1]
# aRA[m]     <- betaRES[[2]]
# aRbetaA[m] <- mean(aRA[1:m])
# #
# betaRES <- betaBMH(nbar = Nbar,
#                    betaAC = betaA, betaBC = betaB,
#                    betaPC = betaP, betaQC = betaQ,
#                    Xa = XA, Xb = XB, Xp = XP, Xq = XQ,
#                    Z = zT, k = kT, VCM = VCMestB)
# betaB      <- betaRES[[1]][1]
# aRB[m]     <- betaRES[[2]]
# aRbetaB[m] <- mean(aRB[1:m])
# #
# betaRES <- betaPMH(nbar = Nbar,
#                    betaAC = betaA, betaBC = betaB,
#                    betaPC = betaP, betaQC = betaQ,
#                    Xa = XA, Xb = XB, Xp = XP, Xq = XQ,
#                    Z = zT, k = kT, VCM = VCMestP)
# betaP      <- betaRES[[1]][1]
# aRP[m]     <- betaRES[[2]]
# aRbetaP[m] <- mean(aRP[1:m])
# #
# betaRES <- betaQMH(nbar = Nbar,
#                    betaAC = betaA, betaBC = betaB,
#                    betaPC = betaP, betaQC = betaQ,
#                    Xa = XA, Xb = XB, Xp = XP, Xq = XQ,
#                    Z = zT, k = kT, VCM = VCMestQ)
# betaQ      <- betaRES[[1]][1]
# aRQ[m]     <- betaRES[[2]]
# aRbetaQ[m] <- mean(aRQ[1:m])