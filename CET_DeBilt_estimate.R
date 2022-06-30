rm(list = ls())
cat("\014")

library(MASS)

################################################################################
############################## Auxiliary Functions #############################
################################################################################

####################################
# Define state space model
####################################

specifyStateSpaceModel = function(hyperParameters, data){
  seasons <- 12
  splineKnots <- c(0, 1, 4, 7, 10, 12)
  k <- length(splineKnots) - 1
  splineWeights <- matrix(0, nrow = seasons, ncol = k)
  nStates <- 2 
  for(i in 1:seasons){
    weights <- periodicSpline(splineKnots, i)
    splineWeights[i, ] <- weights
  }
  eta<-1
  stateSpaceModel <- list(seasons = seasons, splineWeights = splineWeights, data = data, nObs = length(data),
                          eta = eta, seasonalIntercepts = splineWeights %*% hyperParameters[1:k],
                          drifts = splineWeights %*% hyperParameters[(k+1): (2*k)],
                          loadings = splineWeights %*% exp(hyperParameters[(2 * k+1) : (3*k)]), #mapping loadings
                          arCoefs = -1 + 2/(1 + exp(- splineWeights %*% hyperParameters[(3*k+1):(4*k)])), #mapping arCoefs
                          sigmaEps = exp(splineWeights %*% hyperParameters[(4 * k + 1):(5*k)]), #mapping sigmaEps
                          Z = c(0,1), mHH = diag(c(eta,1)), A = numeric(nStates), nStates = nStates,
                          mH=matrix(0,5*k,5*k)) # standard errors
  stateSpaceModel$P <- diag(c(0, stateSpaceModel$sigmaEps[12]/(1 - stateSpaceModel$arCoefs[12]^2)))
  return(stateSpaceModel)
}

####################################
# Estimate state space model
####################################

estimateModel <- function(data,  hyperParameters){
  #switch between Newton-Raphson ("BFGS") and simplex ("Nelder-Mead")
  #optimizer <- optim(hyperParameters, stateSpaceLogLike, data = data, control = list(trace = 4, maxit = 10000, reltol = 1e-9), method="Nelder-Mead", hessian=TRUE)
  optimizer <- optim(hyperParameters, stateSpaceLogLike, data = data, control = list(trace = 4, maxit = 10000, reltol = 1e-9), method="BFGS", hessian=TRUE)
  parameters <- optimizer$par
  stateSpaceModel <- specifyStateSpaceModel(parameters, data)
  stateSpaceModel <- CETKalmanFilter(stateSpaceModel)
  stateSpaceModel <- stateSmoother(stateSpaceModel)
  stateSpaceModel$mu <- mean(stateSpaceModel$seasonalIntercepts)+mean(stateSpaceModel$drifts)*(1:stateSpaceModel$nObs)+mean(stateSpaceModel$loadings)*stateSpaceModel$smoothedA[1,]
  stateSpaceModel$parameters <- parameters
  stateSpaceModel$mH <- optimizer$hessian
  return(stateSpaceModel)
}

########################################
# Evaluate likelihood objective function 
########################################

stateSpaceLogLike <- function(hyperParameters, data){
  stateSpaceModel <- specifyStateSpaceModel(hyperParameters, data)
  seasons <- stateSpaceModel$seasons
  permutationMatrix <- cbind(c(rep(0,seasons-1), 1), diag(seasons)[,-seasons])
  selectionVector <- c(1, numeric(seasons-1))
  intercepts <- stateSpaceModel$seasonalIntercepts
  drifts <- stateSpaceModel$drifts
  loadings <- stateSpaceModel$loadings
  arCoefs <- stateSpaceModel$arCoefs
  sigmaEps <- stateSpaceModel$sigmaEps
  P <- stateSpaceModel$P
  A <- stateSpaceModel$A
  GG <- 0
  llk <- 0
  sumOfSquares <- 0
  for (t in 1:length(data)) {
    Z <- c(t(selectionVector) %*% loadings, 1)
    X <- c(selectionVector, t * selectionVector)
    mT <- diag(c(1, t(selectionVector) %*% arCoefs))
    HH <- diag(c(1, t(selectionVector) %*% sigmaEps))
    if(!is.na(data[t])){ #Observation
      v <- data[t] - Z %*% A - X %*% c(intercepts, drifts)
      f <- Z %*% P %*% Z  + GG
      K <- mT %*% P %*% Z  %*% (1/f) 
      A <- mT %*% A + K %*% v
      P <- mT %*% P %*% mT + HH - (K %*% t(K)) * drop(f)
      llk <- llk + log(f)
      sumOfSquares <- sumOfSquares + v^2/f
    } else { #Missing observation
      A <- mT %*% A
      P <- mT %*% P %*% t(mT) + HH
    }
    selectionVector = t(permutationMatrix) %*% selectionVector
  }
  return((0.5 * (llk + sumOfSquares))[1])
}

####################################
# Run Kalman filter recursions
####################################

CETKalmanFilter <- function(stateSpaceModel){
  data <- stateSpaceModel$data
  nObs <- length(data)
  nStates <- stateSpaceModel$nStates
  seasons <- stateSpaceModel$seasons
  permutationMatrix <- cbind(c(rep(0,seasons-1), 1), diag(seasons)[,-seasons])
  selectionVector <- c(1, numeric(seasons-1))
  llk <- 0
  sumOfSquares <- 0
  stateSpaceModel$innovations = stateSpaceModel$innovationsVar <- rep(NA, nObs)
  stateSpaceModel$statePrediction = stateSpaceModel$statePredictionVar <- matrix(NA, nrow = nObs, ncol = nStates)
  A <- stateSpaceModel$A
  P <- stateSpaceModel$P
  GG <- 0
  for (t in 1:length(data)) {
    Z <- c(t(selectionVector) %*% stateSpaceModel$loadings, 1)
    X <- c(selectionVector, t * selectionVector)
    mT <- diag(c(1, t(selectionVector) %*% stateSpaceModel$arCoefs))
    HH <- diag(c(1, t(selectionVector) %*% stateSpaceModel$sigmaEps))
    if(!is.na(data[t])){ #Observation
      v <- data[t] - Z %*% A - X %*% c(stateSpaceModel$seasonalIntercepts, stateSpaceModel$drifts)
      f <- Z %*% P %*% Z  + GG
      K <- mT %*% P %*% Z  %*% (1/f) 
      A <- mT %*% A + K %*% v
      P <- mT %*% P %*% mT + HH - (K %*% t(K)) * drop(f)
      llk <- llk + log(f)
      sumOfSquares <- sumOfSquares + v^2/f
      stateSpaceModel$innovations[t] <- v
      stateSpaceModel$innovationsVar[t] <- f
    } else { #Missing observation
      A <- mT %*% A
      P <- mT %*% P %*% t(mT) + HH
    }
    selectionVector <- t(permutationMatrix) %*% selectionVector
    stateSpaceModel$statePrediction[t, ] <- A
    stateSpaceModel$statePredictionVar[t, ] <- diag(P)
  }
  stateSpaceModel$llk <- 0.5 * (llk + sumOfSquares)
  return(stateSpaceModel)
}

####################################
# Run state smoothing recursions 
####################################

stateSmoother <- function(stateSpaceModel){
  data <- stateSpaceModel$data
  nObs <- length(data)
  nStates <- stateSpaceModel$nStates
  seasons <- stateSpaceModel$seasons
  permutationMatrix <- cbind(c(rep(0,seasons-1), 1), diag(seasons)[,-seasons])
  selectionVector <- c(1, numeric(seasons-1))
  A <- stateSpaceModel$A
  P <- stateSpaceModel$P
  aV = aF <- array(NA, dim = c(1,1, nObs))
  aaF = aK <- array(NA, dim = c(nStates, 1, nObs))
  aP <- array(NA, dim = c(nStates, nStates, nObs))
  #Kalman recursions
  GG <- 0
  for (t in 1:nObs) {
    Z <- c(t(selectionVector) %*% stateSpaceModel$loadings, 1)
    X <- c(selectionVector, t * selectionVector)
    mT <- diag(c(1, t(selectionVector) %*% stateSpaceModel$arCoefs))
    HH <- diag(c(1, t(selectionVector) %*% stateSpaceModel$sigmaEps))
    aaF[,,t] <- A
    aP[,,t] <- P
    if(!is.na(data[t])){
      v <- data[t] - Z %*% A - X %*% c(stateSpaceModel$seasonalIntercepts, stateSpaceModel$drifts)
      f <- Z %*% P %*% Z  + GG
      K <- mT %*% P %*% Z  %*% (1/f) 
      A <- mT %*% A + K %*% v
      P <- mT %*% P %*% mT + HH - (K %*% t(K)) * drop(f)
      aV[,, t] <- v
      aF[,, t] <- f
      aK[,, t] <- K
    } else { #Missing observation
      A <- mT %*% A
      P <- mT %*% P %*% t(mT) + HH
    }
    selectionVector <- t(permutationMatrix) %*% selectionVector
  }
  #Smoothing recursions
  r <- numeric(nStates)
  n <- matrix(0, nrow = nStates, ncol = nStates)
  stateSpaceModel$smoothedP = stateSpaceModel$smoothedA <- matrix(NA, nrow = nStates, ncol = nObs)
  for (t in nObs : 1) {
    selectionVector <- permutationMatrix %*% selectionVector
    Z <- c(t(selectionVector) %*% stateSpaceModel$loadings, 1)
    mT <- diag(c(1, t(selectionVector) %*% stateSpaceModel$arCoefs))
    if(!is.na(data[t])){
      fInv <- 1/aF[,, t]
      L <- mT - aK[,, t] %*% t(Z)
      r <- Z %*% t(fInv) %*% aV[,, t, drop = FALSE] + t(L) %*% r
      n <- Z %*% t(fInv) %*% Z + t(L) %*% n %*% L
    } else {
      r <- mT %*% r
      n <- t(mT) %*% n %*% mT
    }
    stateSpaceModel$smoothedA[, t] <- aaF[,,t] + aP[,, t] %*% r
    stateSpaceModel$smoothedP[, t] <- diag(aP[,, t] - aP[,, t] %*% n %*% aP[,, t])
  }
  return(stateSpaceModel)
}

####################################
# Calculate periodic spline matrix 
####################################

periodicSpline <- function(x, dX){
  # init
  k <- length(x) - 1
  h <- diff(x)
  h <- c(h, h[1])
  p <- 2 * diag(k)
  q <- matrix(0, nrow = k, ncol = k)
  for (i in 1:(k-1)) {
    p[i, i+1] <- h[i + 1]/(h[i] + h[i + 1])
    p[i+1, i] <- h[i + 1]/(h[i+1] + h[i + 2])
    q[i,i] <- -6/(h[i] * h[i+1])
    q[i, i+1] <- 6 / (h[i+1] * ((h[i] + h[i + 1])))
    q[i+1, i] <- 6/ (h[i+1] * (h[i+1] + h[i + 2]))
  }
  p[1, k] <- h[1] / (h[1] + h[2])
  p[k, 1] <- h[1] / (h[1] + h[k])
  q[1, k] <- 6/(h[1] * (h[1] + h[2]))
  q[k, 1] <- 6 / (h[1] * (h[1] + h[k]))
  q[k, k] <- -6 / (h[1] * h[k])
  if(sum(dX == x) == 1){
    s <- numeric(k)
    if(dX == x[1]){
      r <- c(rep(0, k-1), 1)
    }
    else{
      r <- as.numeric(dX == x[2:length(x)])
    }
  } else {
    ind <- which(dX < x)[1]
    dX_j <- x[ind]
    dX_jMinus1 <- x[ind-1]
    s = r <- numeric(k)
    if (ind == 2){
      r[k] <- (dX_j - dX)/(dX_j - dX_jMinus1)
      r[1] <- (dX - dX_jMinus1)/(dX_j - dX_jMinus1)
      s[k] <- (dX_j - dX) * ( (dX_j - dX)^2 - (dX_j - dX_jMinus1)^2) / (6*(dX_j-dX_jMinus1))
      s[ind - 1] <- (dX - dX_jMinus1) * ( (dX - dX_jMinus1)^2 - (dX_j - dX_jMinus1)^2) / (6 *(dX_j - dX_jMinus1))
    }else if( (ind > 2) && (ind <= k+1)){
      r[ind - 2] <- (dX_j - dX)/(dX_j - dX_jMinus1)
      r[ind - 1] <- (dX- dX_jMinus1)/(dX_j - dX_jMinus1)
      s[ind - 2] <- (dX_j - dX) * ( (dX_j-dX)^2 - (dX_j - dX_jMinus1)^2)/(6*(dX_j - dX_jMinus1))
      s[ind - 1] <- (dX - dX_jMinus1) * ( (dX-dX_jMinus1)^2 - (dX_j - dX_jMinus1)^2)/ (6 * (dX_j - dX_jMinus1))
    }
  }
  weights <- r + s %*% solve(p) %*% q
  return(weights)
}


################################################################################
####################################### MAIN ###################################
################################################################################

########################################
########### CET ########################
########################################

# read data
conn <- url("https://raw.githubusercontent.com/chisos9/data_public/master/cetml1659on_nov_2020.txt")
CET <- read.table(conn, sep = "", skip = 6, header = TRUE,
                  fill = TRUE, na.string = c(-99.99, -99.9))
CET<-CET[96:nrow(CET),] # after 1754
rm(conn)
## store annual mean as separate data frame and remove from data set
Years <- as.numeric(row.names(CET))
ann_mean_CET <- data.frame(MeanTemp = CET[, ncol(CET)],
                           Year = Years)
CET <- CET[, -ncol(CET)]
## "stack" = stack all columns, transforming row names (years) into factor
CET <- stack(CET)[,2:1]
names(CET) <- c("Month","deg C")
## add in Year and nMonth for numeric month and a proper Date class
CET<- transform(CET, Year = (Year <- rep(Years, times = 12)),
                Month = (Month <- rep(1:12, each = length(Years))),
                Date = as.Date(paste(Year, Month, "15", sep = "-"),format = "%Y-%m-%d"))
## sort into temporal order
CET <- CET[with(CET, order(Date)), ]
row.names(CET) <- 1:nrow(CET)
data <- CET$deg.C

# estimate model
vP<-matrix(c(2.6364503666, 7.8897124179, 16.1293937413, 9.3260138651, 3.0110742683, 0.0005924336, 0.0003307760, 0.0001949931, 0.0005922216, 0.0006470174, -4.5701406445, -3.5782102489, -3.7635885898, -3.6156864796, -4.4655601670, 0.6714042120, 0.1362521544, 0.7114851103, 0.1302760087, 0.5434089808, 0.9841380001, 0.2137909968, 0.0741447382, 0.5815697504, 1.3048295542), 25, 1)
stateSpaceModelCET <- estimateModel(data, vP)
vP<-stateSpaceModelCET$parameters
stateSpaceModelCET$Year<-CET$Year

# standard errors by delta-method
mH_<-ginv(stateSpaceModelCET$mH) # Moore-Penrose inverse from MASS package
cS<-stateSpaceModelCET$seasons
mW<-stateSpaceModelCET$splineWeights
ck<-ncol(stateSpaceModelCET$splineWeights)
mJ<-matrix(0,5*cS,5*ck)
mJ[1:cS,1:ck]<-mW
mJ[(cS+1):(2*cS),(ck+1):(2*ck)]<-mW
mJ[(2*cS+1):(3*cS),(2*ck+1):(3*ck)]<-mW%*%exp(vP[(2*ck+1):(3*ck)])
for (i in 1:ck){
  mJ[(3*cS+1):(4*cS),(3*ck+i)]<-2*mW[,i]*exp(-mW%*%vP[(3*ck+1):(4*ck)])/(1+exp(-mW%*%vP[(3*ck+1):(4*ck)]))^2
  mJ[(4*cS+1):(5*cS),(4*ck+i)]<-mW[,i]*exp(mW%*%vP[(4*ck+1):(5*ck)])
}
stateSpaceModelCET$vSe<-sqrt(diag(mJ%*%mH_%*%t(mJ)))

# store results
ResultCET.table<-data.frame(
  mu = stateSpaceModelCET$seasonalIntercepts,
  se_mu = stateSpaceModelCET$vSe[1:12],
  beta = stateSpaceModelCET$drifts,
  se_beta = stateSpaceModelCET$vSe[13:24],
  theta = stateSpaceModelCET$loadings,
  se_theta = stateSpaceModelCET$vSe[25:36],
  phi = stateSpaceModelCET$arCoefs,
  se_phi = stateSpaceModelCET$vSe[37:48],
  nu = stateSpaceModelCET$sigmaEps,
  se_nu = stateSpaceModelCET$vSe[49:60]
)

# calculate phase from Eq. (6) in the paper
seasons = months <- stateSpaceModelCET$seasons #Seasons and months
cK = seasons * months
timestamps <- rep(Years, each = 12) + rep(0:11)/12
sine1 <- sin(2 * pi * 1:seasons/seasons)
cosine1 <- cos(2 * pi * 1:seasons/seasons)
omega <- 2 * pi/cK
c2_CET <- sum(stateSpaceModelCET$seasonalIntercepts * cosine1)/12 + sum(stateSpaceModelCET$drifts * 1:12 * cosine1)/12
c1_CET <- sum(stateSpaceModelCET$seasonalIntercepts * sine1)/12 + sum(stateSpaceModelCET$drifts * 1:12 * sine1)/12
a2_CET <- sum(stateSpaceModelCET$drifts * cosine1)
a1_CET <- sum(stateSpaceModelCET$drifts * sine1)
nYears <- stateSpaceModelCET$nObs/seasons
stateSpaceModelCET$eq6 = atan2(-(c1_CET + a1_CET*(1:(nYears-months))), c2_CET + a2_CET*(1:(nYears-months)))
stateSpaceModelCET$arcSec = (stateSpaceModelCET$eq6[length(stateSpaceModelCET$eq6)] - stateSpaceModelCET$eq6[1])/(nYears - months) * 360 * 60 * 60 /(2*pi)

########################################
########### DeBilt #####################
########################################

# read data
conn <- url("https://raw.githubusercontent.com/chisos9/data_public/master/ilabrijn_nov_2020.txt")
DeBilt1 <- read.table(conn, sep = "", skip = 15, header = FALSE,
                      fill = TRUE, na.string = c("-999.9000"))
DeBilt<-DeBilt1[,-1]
rownames(DeBilt)<-DeBilt1[,1]
colnames(DeBilt)<-c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
rm(conn)
## store annual mean as separate data frame and remove from data set
Years <- as.numeric(row.names(DeBilt))
ann_mean_DeBilt <- data.frame(MeanTemp = DeBilt[, ncol(DeBilt)],
                              Year = Years)
## "stack" = stack all columns, transforming row names (years) into factor
DeBilt <- stack(DeBilt)[,2:1]
names(DeBilt) <- c("Month","deg C")
## add in Year and nMonth for numeric month and a proper Date class
DeBilt<- transform(DeBilt, Year = (Year <- rep(Years, times = 12)),
                   Month = (Month <- rep(1:12, each = length(Years))),
                   Date = as.Date(paste(Year, Month, "15", sep = "-"),format = "%Y-%m-%d"))
## sort into temporal order
DeBilt <- DeBilt[with(DeBilt, order(Date)), ]
row.names(DeBilt) <- 1:nrow(DeBilt)
data <- DeBilt$deg.C

# estimate model
vP<-matrix(c(1.0905864929,8.4547675738,17.2645546569,10.6002709281,2.1740963728,0.0006709707,0.0004305480,0.0003204317,0.0002656379,0.0006695093,-3.3643526312,-3.2768097918,-3.5292404294,-3.3731564271,-3.1786987285,0.6999287763,0.1253102833,0.5585520713,0.2665742531,0.6193863136,1.4749460744,0.6181701467,0.3056625219,0.9087482876,1.9467557645),25,1)
stateSpaceModelDeBilt <- estimateModel(data, vP)
vP<-stateSpaceModelDeBilt$parameters
stateSpaceModelDeBilt$Year<-DeBilt$Year

# standard errors by delta-method
mH_<-ginv(stateSpaceModelDeBilt$mH) # Moore-Penrose inverse from MASS package
cS<-stateSpaceModelDeBilt$seasons
mW<-stateSpaceModelDeBilt$splineWeights
ck<-ncol(stateSpaceModelDeBilt$splineWeights)
mJ<-matrix(0,5*cS,5*ck)
mJ[1:cS,1:ck]<-mW
mJ[(cS+1):(2*cS),(ck+1):(2*ck)]<-mW
mJ[(2*cS+1):(3*cS),(2*ck+1):(3*ck)]<-mW%*%exp(vP[(2*ck+1):(3*ck)])
for (i in 1:ck){
  mJ[(3*cS+1):(4*cS),(3*ck+i)]<-2*mW[,i]*exp(-mW%*%vP[(3*ck+1):(4*ck)])/(1+exp(-mW%*%vP[(3*ck+1):(4*ck)]))^2
  mJ[(4*cS+1):(5*cS),(4*ck+i)]<-mW[,i]*exp(mW%*%vP[(4*ck+1):(5*ck)])
}
stateSpaceModelDeBilt$vSe<-sqrt(diag(mJ%*%mH_%*%t(mJ)))

# store results
ResultDeBilt.table<-data.frame(
  mu = stateSpaceModelDeBilt$seasonalIntercepts,
  se_mu = stateSpaceModelDeBilt$vSe[1:12],
  beta = stateSpaceModelDeBilt$drifts,
  se_beta = stateSpaceModelDeBilt$vSe[13:24],
  theta = stateSpaceModelDeBilt$loadings,
  se_theta = stateSpaceModelDeBilt$vSe[25:36],
  phi = stateSpaceModelDeBilt$arCoefs,
  se_phi = stateSpaceModelDeBilt$vSe[37:48],
  nu = stateSpaceModelDeBilt$sigmaEps,
  se_nu = stateSpaceModelDeBilt$vSe[49:60]
)

# calculate phase from Eq. (6) in the paper
seasons = months <- stateSpaceModelDeBilt$seasons #Seasons and months
cK = seasons * months
timestamps <- rep(Years, each = 12) + rep(0:11)/12
sine1 <- sin(2 * pi * 1:seasons/seasons)
cosine1 <- cos(2 * pi * 1:seasons/seasons)
omega <- 2 * pi/cK
c2_CET <- sum(stateSpaceModelDeBilt$seasonalIntercepts * cosine1)/12 + sum(stateSpaceModelDeBilt$drifts * 1:12 * cosine1)/12
c1_CET <- sum(stateSpaceModelDeBilt$seasonalIntercepts * sine1)/12 + sum(stateSpaceModelDeBilt$drifts * 1:12 * sine1)/12
a2_CET <- sum(stateSpaceModelDeBilt$drifts * cosine1)
a1_CET <- sum(stateSpaceModelDeBilt$drifts * sine1)
nYears <- stateSpaceModelDeBilt$nObs/seasons
stateSpaceModelDeBilt$eq6 = atan2(-(c1_CET + a1_CET*(1:(nYears-months))), c2_CET + a2_CET*(1:(nYears-months)))
stateSpaceModelDeBilt$arcSec = (stateSpaceModelDeBilt$eq6[length(stateSpaceModelDeBilt$eq6)] - stateSpaceModelDeBilt$eq6[1])/(nYears - months) * 360 * 60 * 60 /(2*pi)