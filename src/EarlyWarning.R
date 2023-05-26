
#########################################################################################
# Function to extract forecast errors and forecast variances from a result object of a 
# DLM function. Some times they are stored as lists and sometimes they are stored as vectors.
# The function assumes that the errors are univariate.
#########################################################################################
# res: The result object returned by a DLM function
#########################################################################################
# Returna a list with two elements:
# - a vector of forecast errors
# - a vector of forecast variances
#########################################################################################
extractErrorsFromResult = function(res) {
  errors = c()
  variances = c()
  # Are the forecast errors stored as a list or a vector
  if (is.list(res$et)) {
    # They are stored as a list
    errors = unlist(res$et)
    variances = unlist(res$Qt)
  }
  else {
    # They are stored as a vector
    errors = res$et
    variances = res$Qt
  }
  return(list(et = errors, Qt = variances))
}


#########################################################################################
# Function to create a Shewhart control chart for univariate models.
#########################################################################################
# res: The result object returned by a DLM function
# limits: The number of units of standard deviations for the control limits
#########################################################################################
runShewhartUnivariate = function(res, limits) {
  errors = extractErrorsFromResult(res)
  minimum = 1.1*min(errors$et, na.rm = TRUE)
  maximum = 1.1* max(errors$et, na.rm = TRUE)
  plot(errors$et, type = "l", lwd= 2, ylim = c(minimum, maximum))
  abline(h=0)
  lines(limits*sqrt(errors$Qt))
  lines(-limits*sqrt(errors$Qt))
}




#########################################################################################
# Function to calculate the cumulative sum of forecast errors for univariate models.
# The forecast errors can, optionally, be standardized.
#########################################################################################
# res: The result object returned by a DLM function
# standardize: Must be TRUE if the errors should be standardized
#########################################################################################
runCusumUnivariate = function(res, standardize = FALSE) {
  errors = extractErrorsFromResult(res)
  if (standardize) {
    errors$et = errors$et/sqrt(errors$Qt)
  }
  cusum = c()
  cusum[1] = 0
  for (i in 2:length(errors$et)) {
    cusum[i] = sum(errors$et[2:i])
  }
  return(cusum)
}


##############################################################################
# Runs a V-mask on the Cusum. The function is later used by the function     #
# runAndPlotV.mask. Refer also to Figure 8.28 and the description on Pages   #
# 166-167.                                                                   #
##############################################################################
# Use res as input. It is output from the runDGLM function.                  #
# Use dist as the lead distance                                              #
# Use angle as the slope of the arms of the V-mask.                          #
# Use reset as a flag for resetting the Cusum to zero after an alarm.        #
##############################################################################
runV.mask = function(res, dist, angle, reset) {
  errors = extractErrorsFromResult(res)
  cusums = runCusumUnivariate(res, TRUE)
  mask = matrix(rep(F, length(errors$et)*6), length(errors$et))
  mask[1, 1] = 1 
  mask[1, 2] = cusums[1]
  mask[1, 6] = cusums[1]
  lastAlarm = 0
  for (i in 2:length(errors$et)) {
    mask[i, 1] = i 
    mask[i, 2] = cusums[i]
    alarmLastWeek = (mask[i-1, 3] || mask[i-1, 4])
    # Reset to 0 if alarm last week
    if (alarmLastWeek) {
      mask[i, 6] = mask[i, 2] - mask[i-1, 2]
    }
    else {
      mask[i, 6] = mask[i, 2] - mask[i-1, 2] + mask[i-1, 6]        
    }
    currentLevel =  cusums[i]
    if (reset) currentLevel = mask[i, 6]
    if (reset == F || alarmLastWeek == F) {
      for (j in (i-1):(lastAlarm+1)) {
        point = i + dist
        lowerArm = currentLevel - (i-j+dist)*angle
        upperArm = currentLevel + (i-j+dist)*angle
        previousLevel = cusums[j]
        if (reset) previousLevel = mask[j, 6]      
        if (previousLevel < lowerArm) {
          mask[i, 3] = T
          mask[i, 5] = j
          if (reset) lastAlarm = i
        }
        if (previousLevel > upperArm) {
          mask[i, 4] = T
          mask[i, 5] = j
          if (reset) lastAlarm = i
        }
      }
    }
  }
  return(mask)
}


##############################################################################
# Adds a V-mask to an existing Cusum plot. The function is later used by the #
# function runAndPlotV.mask.                                                 #
##############################################################################
# Use mask as input. It is output from the runV.mask function.               #
# Use week as the week to plot the V-mask                                    #
# Use dist as the lead distance                                              #
# Use angle as the slope of the arms of the V-mask.                          #
# Use reset as a flag for resetting the Cusum to zero after an alarm.        #
##############################################################################
addMaskToPlot = function(mask, week, dist, angle, reset) {
  size = min(20, week-1)
  rows = max(length(mask[,1]), week+dist)
  arms = matrix(rep(NA, rows*4), rows)
  base = mask[week, 2]
  if (reset) base = mask[week, 6]
  for (i in (week+dist):(week+dist-size)) {
    arms[i,1] = base - (week + dist - i) * angle 
    if (i > week - 1) arms[i, 2] = base
    arms[i,3] = base + (week + dist - i) * angle 
    arms[i, 4] = i
  }  
  lines(arms[,4], arms[,1], lty = 2, col = "red", lwd = 1)
  lines(arms[,4], arms[,2], lty = 2, col = "red", lwd = 1)
  lines(arms[,4], arms[,3], lty = 2, col = "red", lwd = 1)
}

##############################################################################
# Runs a V-mask on the Cusum, identifies all alarms, plots the Cusum with a  #
# V-mask for each alarm. Refer also to Figure 8.28 and the description on    #
# Pages 166-167. The function first calls the runV.mask function, then it    #
# plots the Cusum, identifies all alarms and adds V-masks to the plot.       #
##############################################################################
# Use res as input. It is output from the runDGLM function.                  #
# Use dist as the lead distance                                              #
# Use angle as the slope of the arms of the V-mask.                          #
# Use ylims to define the scope of the y axis                                #
# Use reset as a flag for resetting the Cusum to zero after an alarm.        #
##############################################################################
runAndPlotV.mask = function(res, dist, angle, ylims, reset) {
  mask = runV.mask(res, dist, angle, reset)
  # The Cusum plot
  if (reset) {
    plot(mask[,1], mask[,6], t = 'l', xlab = "Time", ylab="Cumulated sum of standardized forecast errors", col = "blue", ylim = ylims)
  }
  else {
    plot(mask[,1], mask[,2], t = 'l', xlab = "Time", ylab="Cumulated sum of standardized forecast errors", col = "blue", ylim = ylims)  
  }
  abline(h = 0.0)
  for (i in 1:length(mask[,1])) {
    if (mask[i, 3] || mask[i, 4]) addMaskToPlot(mask, i, dist, angle, reset)
  }
  return(mask)
}



############################################################################################
# Function calculating the (upper and lower) tabular cusum for forecast errors
############################################################################################
# res: The output object of a a runDLM function
# slack: The slack value
# standardized: Use TRUE if the forecast errors should be standardized
############################################################################################
# Returns a list with two vectors of cusums (upper and lower)
############################################################################################
runTabularCusumUnivariate = function(res, slack, standardize = FALSE) {
  errors = extractErrorsFromResult(res)
  if (standardize) {
    errors$et = errors$et/sqrt(errors$Qt)
  }
  upperC = c()
  lowerC = c()
  upperC[1] = 0
  lowerC[1] = 0
  for (i in 2:length(errors$et)) {
    upperC[i] = max(0, errors$et[i]- slack + upperC[i-1])
    lowerC[i] = max(0, -slack-errors$et[i] + lowerC[i-1])
  }
  return(list(upper = upperC, lower = lowerC))
}

############################################################################################
# Function creating a cusum plot with decision lines
############################################################################################
# res: The output object of a a runDLM function
# slack: The slack value
# decision: The decision level
# standardized: Use TRUE if the forecast errors should be standardized
############################################################################################
createTabularCusumPlot = function(res, slack, decision, standardize = FALSE) {
  tabs = runTabularCusumUnivariate(dlm, slack, TRUE)
  amp = max(abs(tabs$upper), abs(tabs$lower))
  plot(tabs$upper, type ='h', ylim = c(-amp, amp))
  lines(-tabs$lower, type='h')
  abline(h=decision, col = "red")
  abline(h=-decision, col = "red")
  abline(h = 0)
}


############################################################################################
# Function returning Mahalanobis distance
############################################################################################
# mean: The mean vector
# variance: The variance-covariance matrix
############################################################################################
# Returns Mahalanobis distance
############################################################################################
findMahalanobisDistanceFromZero = function(mean, variance) {
  return(sqrt(t(mean) %*% solve(variance) %*% mean))
}


#########################################################################################
# Function to extract Mahalanobis distances from a result object of a  multivariate DLM 
# The mean and variances of the forecasts are stored as lists.
#########################################################################################
# res: The result object returned by a multivariate DLM function
#########################################################################################
# Returns a vector of Mahalanobis distances from zero
#########################################################################################
extractMahalanobisDistances = function(res)  {
  mDists = c()
  n = length(res$et)
  for (i in 1:n) {
    mDists[i] = findMahalanobisDistanceFromZero(res$et[[i]], res$Qt[[i]])
  }
  return(mDists)
}

#########################################################################################
# Function to extract squares of Mahalanobis distances from a result object of a  
# multivariate DLM The mean and variances of the forecasts are stored as lists.
#########################################################################################
# res: The result object returned by a multivariate DLM function
#########################################################################################
# Returns a vector of squares of Mahalanobis distances from zero
#########################################################################################
extractSquaredMahalanobisDistances = function(res)  {
  mDists = c()
  n = length(res$et)
  for (i in 1:n) {
    mDists[i] = findMahalanobisDistanceFromZero(res$et[[i]], res$Qt[[i]])^2
  }
  return(mDists)
}

#########################################################################################
# Function to create a generalized Shewhart control chart for multivariate models.
# The squares of Mahalinobis distances is plottet together with an upper control limit
# for the chi-square distribution.
#########################################################################################
# res: The result object returned by a multivariate DLM function
# quantile: The quantile of the upper control limit
#########################################################################################
runChiSquarePlot = function(res, quantile) {
  squares = extractSquaredMahalanobisDistances(res)
  freedom = length(res$et[[1]])
  plot(squares, type ='l')
  abline(h=qchisq(quantile, freedom))
}
