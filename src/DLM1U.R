##############################################################
# FUNCTIONS FOR FIRST ORDER UNIVARIATE POLYNOMIAL DLMs			 #
##############################################################

##############################################################
# Functiuon running a first order univariate DLM
##############################################################
# Yt: The data series to be monitored
# mu0: Prior mean
# C0: Prior variance
# V: Observation variance
# W: System variance
##############################################################
# Returns results as a list
##############################################################
DLM1U <- function(Yt,mu0,C0,V,W){
  
  n <- length(Yt);
  at <- c(); 			        # Define the vectors
  Rt <- c();
  ft <- c();
  Qt <- c();
  At <- c();
  et <- c();
  mt <- c();
  Ct <- c();
  
  mt[1] <- mu0;				# Prior Distribution
  Ct[1] <- C0;
  
  for(i in (2:n)) {
    
    at[i] <- mt[i-1];		    		# Prior mean
    Rt[i] <- Ct[i-1] + W ;			# Prior Variance
    
    ft[i] <- at[i];				# One-step Forecast mean
    Qt[i] <- Rt[i] + V;			# One-step Forecast variance
    
    At[i] <- Rt[i] / Qt[i];			# Adaptative Coef. matrix
    et[i] <- Yt[i] - ft[i];			# one-step forecast error
    
    # Check for missing data
    if (!is.na(Yt[i])) {
      # Not missing - update as described in algorithm
      mt[i] <- at[i] + At[i]*et[i];		# Filtered mean
    }
    else {
      # Missing: Do not update estimated mean
      mt[i] = mt[i-1]
      et[i] = 0
    }
    Ct[i] <- At[i] * V;			# Filtered variance - update even if missing
    
  }
  
  return(list(
    at=at,
    Rt=Rt,
    ft=ft,
    Qt=Qt,
    At=At,
    et=et,
    mt=mt,
    Ct=Ct))
}


##############################################################
# Functiuon running a first order univariate DLM with discounting
##############################################################
# Yt: The data series to be monitored
# mu0: Prior mean
# C0: Prior variance
# V: Observation variance
# delta: The discount factor
##############################################################
# Returns results as a list
##############################################################
DLM1UD = function(Yt,mu0,C0,V,delta){

	n <- length(Yt);
	at <- c(); 			        # Define the vectors
	Rt <- c();
	ft <- c();
	Qt <- c();
	At <- c();
	et <- c();
	mt <- c();
	Ct <- c();
	
	mt[1] <- mu0;				# Prior Distribution
 	Ct[1] <- C0;
 	et[1] = 0
	
  for(i in (2:n)) {
	
	  at[i] <- mt[i-1];		    	# Prior mean
	  Rt[i] <- Ct[i-1] / delta ;		# Prior Variance
	
	  ft[i] <- at[i];				# One-step Forecast mean
	  Qt[i] <- Rt[i] + V;			# One-step Forecast variance
 
	  At[i] <- Rt[i] / Qt[i];			# Adaptative Coef. matrix
	  et[i] <- Yt[i] - ft[i];			# one-step forecast error
	
	  if (!is.na(Yt[i])) {
	    mt[i] <- at[i] + At[i]*et[i];		# Filtered mean
	  }
	  else {
  	  mt[i] = mt[i-1]
	    et[i] = 0
	  }
	  Ct[i] <- At[i] * V;			# Filtered variance
	
  }

	return(list(
		at=at,
		Rt=Rt,
		ft=ft,
		Qt=Qt,
		At=At,
		et=et,
		mt=mt,
		Ct=Ct));
	};


################################################################################
# Function running the Kalman smoother for a first order univariate DLM 
################################################################################
# res: The result of either DLM1U or DLM1UD
################################################################################
# Returns results as a matrix - first column means, second variances
################################################################################
smoother1U = function(res) {
  
  n = length(res$mt)

  sm = matrix(NA, n, 2)
  
  # Put last value equal to filtered
  sm[n, 1] <- res$mt[n]
  sm[n, 2] <- res$Ct[n]
  
  # Iterate backwards over days
  for(i in ((n-1):1))   {
    
    
    Bt <- res$Ct[i] / res$R[i+1]
    sm[i, 1] <- res$mt[i] + Bt * (sm[i+1,1] - res$a[i+1])
    sm[i, 2] <- res$Ct[i] + Bt * (sm[i+1,2] - res$R[i+1]) * Bt
  }
  
  return(sm)
}


################################################################################
# Function running the Kalman filter and returning the sum of squares of forecast errors
################################################################################
# pars: A vecor with two elements (V and W)
# Yt: The data series
# m0: The prior mean
# C0: The prior variance
################################################################################
# Returns the sum of squares of forecast errors
################################################################################
SSSFE = function(pars, Yt, m0, C0) {
  print(paste("Now testing V =", pars[1], "W =", pars[2]))
  res = DLM1U(Yt, m0, C0, pars[1], pars[2])
  sse = (res$et)^2
  print(paste("Now testing V =", pars[1], "W =", pars[2], "SSFE =", sum(sse[2:length(sse)])))
  return(sum(sse[2:length(sse)]))
}


################################################################################
# Function estimating the variance components in a first order univariate 
# polynomial DLM by numerical minimization of the sum of squares of forecast
# errors.
################################################################################
# Yt: The dataseries
# m0: The prior mean
# C0: The prior variance
# startV: The initial guess for the observation variance
# startW: The initial guess for the system variance
################################################################################
# Returns a vector with the optimal variances (V and W)
################################################################################
estimateByOptimization = function(Yt, m0, C0, startV, startW) {
  res = optim(c(startV, startW), SSSFE, Yt=Yt, m0=m0, C0=C0)
  return(res)
}

################################################################################
# Function estimating the variance components of a first order univariate
# polynomial DLM. The EM algorithm as described by Dethlefsen (2001)
################################################################################
# Yt: The data series
# m0: The prior mean
# C0: The prior variance
# startV: The initial guess for the observation variance
# startW: The initial guess for the system variance
# toleranceV: Observation variances of two iterations are considered equal if 
#             they deviate less than this value
# toleranceW: System variances of two iterations are considered equal if 
#             they deviate less than this value
################################################################################
# Returns a vector with the optimal variances (V and W)
################################################################################
runEM1U = function(Yt, m0, C0, startV, startW, 
                   toleranceV = 0.00000001, toleranceW  = 0.00000001) {
  newV = startV
  newW = startW
  oldV = 0
  oldW = 0
  n = length(Yt)
  c = 0
  print(newV)
  print(oldV)
  print(newW)
  print(oldW)
  while(abs(newV - oldV) > toleranceV | abs(newW - oldW) > toleranceW) {
    oldV = newV
    oldW = newW
    res = DLM1U(Yt, m0, C0, oldV, oldW)
    smo = smoother1U(res)
    sumV = 0
    sumW = 0
    c = c + 1
    nonMis = 0
    for (t in 2:(n)) {
      # For V, only include non-missing observations (and count how many included)
      if (!is.na(Yt[t])) {
        sumV = sumV + smo[t, 2] + (Yt[t]-smo[t,1])^2
        nonMis = nonMis +1
      }
      Bt_1 = res$Ct[t-1]/res$Rt[t]
      Lt = smo[t, 2]+ smo[t-1, 2]-2*smo[t, 2]*Bt_1
      sumW = sumW + Lt + (smo[t, 1] - smo[t-1,1])^2
    }
    newV = sumV/nonMis
    newW = sumW/(n-1)
    print(paste("Iteration", c, "with", "V =",newV, "and W =", newW  ))
  }
  return(c(newV, newW))
}


# ALARM FOR CUSUM

FUNALARM2 <- function (cu,lim) {
  n <- length(cu);
  i <- 2;
  a <- 1;
  res <- c();
  while(i < n){
  if(cu[i] >= lim & cu[i-1] < lim) 
	{res[a] <- i; i <- i+1; a <- a+1} else {i <- i+1}
  }
  cat("alarm at hour ", res);
  }


DLMExtended = function(Yt, m0, C0, Ft, Gt, V, W) {
  
  n <- length(Yt);
  at <- list(); 			        # Define the vectors
  Rt <- list();
  ft <- c();
  Qt <- c();
  At <- list();
  et <- c();
  mt <- list();
  Ct <- list();
  
  mt[[1]] <- m0;				# Prior Distribution
  Ct[[1]] <- C0;
  
  for(i in (2:n)) {
    
    at[[i]] <- Gt %*% mt[[i-1]];		    		# Prior mean
    Rt[[i]] <- Gt %*% Ct[[i-1]] + W ;			# Prior Variance
    
    ft[i] <- t(Ft) %*% at[[i]];				# One-step Forecast mean
    Qt[i] <- t(Ft)%*% Rt[[i]] %*% Ft + V;			# One-step Forecast variance
    
    At[[i]] <- Rt[[i]] %*% Ft %*% solve(Qt[i]);			# Adaptative Coef. matrix
    et[i] <- Yt[i] - ft[i];			# one-step forecast error
    
    # Check for missing data
    if (!is.na(Yt[i])) {
      # Not missing - update as described in algorithm
      mt[[i]] <- at[[i]] + At[[i]]*et[i];		# Filtered mean
    }
    else {
      # Missing: Do not update estimated mean
      mt[[i]] = mt[[i-1]]
      et[i] = 0
    }
    Ct[[i]] <- Rt[[i]] - At[[i]]%*%Qt[i]%*% t(At[[i]]);			# Filtered variance - update even if missing
    
  }
  
  return(list(
    at=at,
    Rt=Rt,
    ft=ft,
    Qt=Qt,
    At=At,
    et=et,
    mt=mt,
    Ct=Ct))
}

smootherExtended = function(res, Ft, Gt) {
  
  
  n = length(res$mt)
  sm = matrix(NA, n, 2)
  
  # Put last value equal to filtered
  sm[n, 1] <- sum(res$mt[[n]])
  sm[n, 2] <- t(Ft)%*%res$Ct[[n]]%*%Ft
  oldSm = res$Ct[[n]]
  oldM = res$mt[[n]]
  
  
  # Iterate backwards over days
  for(i in ((n-1):1))   {
    
    
    Bt <- res$Ct[[i]]%*%t(Gt)%*% solve(res$R[[i+1]])
    newM <- res$mt[[i]] + Bt %*% (oldM - res$a[[i+1]])
    newSm <- res$Ct[[i]] + Bt %*% (oldSm - res$R[[i+1]]) %*% t(Bt)
    sm[i, 1] = t(Ft) %*% newM
    sm[i, 2] = t(Ft) %*% newSm %*% Ft
    oldSm = newSm
    oldM = newM
  }
  
  return(sm)
}

  
  
runForecast = function(mt, Ct, V, W, Ft, Gt, steps) {
  print(mt)
  at = list()
  Rt = list()
  ft = list()
  Qt = list()
  
  at[[1]] = Gt %*% mt
  Rt[[1]] = Gt %*% Ct %*% t(Gt) + W
  ft[[1]] = t(Ft) %*% at[[1]]
  Qt[[1]] = t(Ft) %*% Rt[[1]] %*% Ft + V
  
  for (i in 2:steps) {
    at[[i]] = Gt %*% at[[i-1]]
    Rt[[i]] = Gt %*% Rt[[i-1]] %*% t(Gt) + W
    ft[[i]] = t(Ft) %*% at[[i]]
    Qt[[i]] = t(Ft) %*% Rt[[i]] %*% Ft + V
  }
  return(list(at=at, Rt=Rt, ft=ft, Qt=Qt))
}
