#### Posterior Simulator Functions ####
#
#
#
#
#
### Gibbs block for variance (sigma_mu)^(-2) of individual heterogeneity mu:
tauGibbs <- function(a1, s1) {
  s <- 1/(s1)
  sqrt(1/rgamma(n = 1, shape = a1, scale = s))
}
#
#
#
#
#
### MH step for phi
betaMHall <- function(nbars, betaAC, betaBC, betaPC, betaQC,
                      ceptAC, ceptBC, ceptPC, ceptQC,
                      Xa, Xb, Xp, Xq,
                      Z, k, VCM, N, t,
                      indba, indbb, indbp, indbq,
                      indca, indcb, indcp, indcq) {
  accept <- 0
  betaC  <- c(betaAC, betaBC, betaPC, betaQC, ceptAC, ceptBC, ceptPC, ceptQC)
  ##############################################################################
  betaP <- mvrnorm(n = 1, mu = betaC, Sigma = VCM)
  # while (any(betaP <= 0)) {
  #   betaP <-  mvrnorm(n = 1, mu = betaC, Sigma = VCM)
  # }
  ##############################################################################
  AC <- exp(Xa %*% betaAC + ceptAC)
  BC <- exp(Xb %*% betaBC + ceptBC)
  PC <- exp(Xp %*% betaPC + ceptPC)
  QC <- exp(Xq %*% betaQC + ceptQC)
  
  AP <- exp(Xa %*% betaP[indba] + betaP[indca])
  BP <- exp(Xb %*% betaP[indbb] + betaP[indcb])
  PP <- exp(Xp %*% betaP[indbp] + betaP[indcp])
  QP <- exp(Xq %*% betaP[indbq] + betaP[indcq])
  ##############################################################################
  Cpart <- 0
  Ppart <- 0
  for (j in 1:t){
    for (i in 1:N) {
      ind <- (j - 1)*Ncross + i
      zi     <- Z[[ind]][2:k]
      zi1    <- Z[[ind]][1:k]
      ############################################################################
      ui    <- ((zi/BC[ind])^AC[ind])/(1 + (zi/BC[ind])^AC[ind])
      ui1   <- ((zi1/BC[ind])^AC[ind])/(1 + (zi1/BC[ind])^AC[ind])
      Fi  <- pbeta(q = ui, shape1 = PC[ind], shape2 = QC[ind])
      Fi  <- c(Fi, 1)
      Fi1 <- pbeta(q = ui1, shape1 = PC[ind], shape2 = QC[ind])
      
      piC <- Fi - Fi1
      piC <- log(piC)
      
      piC    <- piC*nbars[[ind]]
      Cpart  <- Cpart + sum(piC)
      ##########################################################################
      ui    <- ((zi/BP[ind])^AP[ind])/(1 + (zi/BP[ind])^AP[ind])
      ui1   <- ((zi1/BP[ind])^AP[ind])/(1 + (zi1/BP[ind])^AP[ind])
      Fi  <- pbeta(q = ui, shape1 = PP[ind], shape2 = QP[ind])
      Fi  <- c(Fi, 1)
      Fi1 <- pbeta(q = ui1, shape1 = PP[ind], shape2 = QP[ind])
      
      piP <- Fi - Fi1
      piP <- log(piP)
      
      piP    <- piP*nbars[[ind]]
      Ppart  <- Ppart + sum(piP)
      ##########################################################################
    }
  }
  ##############################################################################
  alpha <- (Ppart - Cpart) 
  if (log(runif(1)) <= alpha) {
    beta   <- betaP
    accept <- 1
  } else {
    beta   <- betaC
    accept <- 0
  }
  ##############################################################################
  return(list(beta, accept))
}
DGP2 <- function(N, A, B, P, Q, Z) {
  y <- rgb2(n = N, shape1 = A, scale = B, shape2 = P, shape3 = Q)
  ncut <- cut(y, breaks = Z)
  ncount <- table(ncut)
  return(as.vector(ncount))
}
lpost <- function(nbars, pars, Xa, Xb, Xp, Xq, Z, k, N, t,
                  indba, indbb, indbp, indbq,
                  indca, indcb, indcp, indcq) {
  A   <- exp(Xa %*% pars[indba] + pars[indca])
  B   <- exp(Xb %*% pars[indbb] + pars[indcb])
  P   <- exp(Xp %*% pars[indbp] + pars[indcp])
  Q   <- exp(Xq %*% pars[indbq] + pars[indcq])
  out <- 0
  for (j in 1:t){
    for (i in 1:N) {
      ind <- (j - 1)*Ncross + i
      zi  <- Z[[ind]][2:k]
      zi1 <- Z[[ind]][1:k]
      ui  <- ((zi/B[ind])^A[ind])/(1 + (zi/B[ind])^A[ind])
      ui1 <- ((zi1/B[ind])^A[ind])/(1 + (zi1/B[ind])^A[ind])
      Fi  <- pbeta(q = ui, shape1 = P[ind], shape2 = Q[ind])
      Fi  <- c(Fi, 1)
      Fi1 <- pbeta(q = ui1, shape1 = P[ind], shape2 = Q[ind])
      
      PreOut <- Fi - Fi1
      PreOut <- log(PreOut)
      PreOut <- PreOut*nbars[[ind]]
      out    <- out + sum(PreOut)
    }
  }
  return(out)
}


# 
# ll <- function(nbars, pars, Z, k) {
#   A <- pars[1]
#   B <- pars[2]
#   P <- pars[3]
#   Q <- pars[4]
#   zi  <- Z[2:k]
#   zi1 <- Z[1:k]
#   ui  <- ((zi/B)^A)/(1 + (zi/B)^A)
#   ui1 <- ((zi1/B)^A)/(1 + (zi1/B)^A)
#   Fi  <- pbeta(q = ui, shape1 = P, shape2 = Q)
#   Fi  <- c(Fi, 1)
#   Fi1 <- pbeta(q = ui1, shape1 = P, shape2 = Q)
#   
#   out <- Fi - Fi1
#   out <- log(out)
#   out <- out*nbars[[i]]
#   out <- sum(out)
#   return(out)
# }