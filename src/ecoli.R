library(dplyr)
library(tidyr)
library(readxl)
library(lubridate)

d <- read_excel("./data/raw/Necropsy.xlsx")

nec.dat <- d %>%
    select(RegistrationDate,
           Sequence,
           DiagnosisCodesFlat,
           AnimalType) %>%
    filter(AnimalType %in% c('SS','SSB','SSR','SSS','SST','SSU','SSV'),
           Sequence == 0) %>%
    mutate(M = month(RegistrationDate),
           Y = (as.numeric(year(RegistrationDate)) - 2013) * 12,
           t = Y + M,
           ECO = case_when(grepl('ECO', 
                                 DiagnosisCodesFlat, 
                                 fixed = T) ~ 1,
                     TRUE ~ 0))
View(nec.dat)

# dECO has:
# 1) p = aggregated prevalence by month (no counts)
# 2) count = number of positive results (= has mention of E. coli)
# 3) n = number of total diagnoses 

dECO <- nec.dat %>%
    group_by(t) %>%
    summarise(p = mean(ECO),
              count = sum(ECO),
              n = n()) %>%
    mutate(log_p = log(p/(1 - p)))
View(dECO)

###############################################

# avoid scientific notification
options(scipen=999)

# use formulas
source('DLM1U.R')

# plots
p_plot <- plot(dECO$p ~ dECO$t, 
     xlab='Year-month', 
     ylab='Fraction E.coli diagnosis')
p_plot

c_plot <- plot(dECO$count ~ dECO$t, 
     xlab='Year-month', 
     ylab='Total E.coli diagnosis')
c_plot
abline(v=c(seq(0.5, 120, 12)), col= "black")


# first order polynomial
# learning set = first 10 months
learn_set <- dECO[2:20,]

# use smoother to get best estimates
mean_count = mean(learn_set$count)
var_count = var(learn_set$count)
mean_p = mean(learn_set$p)
var_p = var(learn_set$p)

# plot counts
dlm = DLM1U(dECO$count, mean_count, var_count, 1, 0.04)
sm1U = smoother1U(dlm)

plot(dECO$count ~ dECO$t, 
     xlab ='Year-month', 
     ylab ='Total E.coli diagnosis',
     type = "l")

lines(dlm$mt, col = "green", lwd = 2)
lines(sm1U[,1], col = "blue", lwd = 2 )

# plot prevalence
dlm_p = DLM1U(dECO$p, mean_p, var_p, 1, 0.04)
sm1U_p = smoother1U(dlm_p)

plot(dECO$p ~ dECO$t, 
     xlab='Year-month', 
     ylab='Prevalence E.coli diagnosis',
     type = 'l')

lines(dlm_p$mt, col = "green", lwd = 2)
lines(sm1U_p[,1], col = "blue", lwd = 2 )

#### improve the function ####

# starting points from smoother function:
mu0 = sm1U_p[1,1]
C0 = sm1U_p[1,2]
    
# Find optimal variance components by EM algorithm
optiEM = runEM1U(learn_set$p, mu0, C0, 1, 0.04)
bestV = optiEM[1]
bestW = optiEM[2]
optiRatioEM = bestV/bestW

# run filter and smoother with new values
dlm = DLM1U(learn_set$p, mu0, C0, bestV, bestW)
sm1U = smoother1U(dlm)
plot(learn_set$p, type="l", ylim = c(0,1),
     ylab = "Prevalence",
     xlab = "Year-month")
lines(dlm$mt, col = "green", lwd = 2)
lines(sm1U[,1], col = "blue", lwd = 2 )
lines((sm1U[,1] + 2*sqrt(sm1U[,2])), col = "blue", lty=2)
lines((sm1U[,1] - 2*sqrt(sm1U[,2])), col = "blue", lty=2)


### assessment ###    

# forecast errors should be normally distributed
n = length(dlm$et)
plot(dlm$et[1:(n-1)], dlm$et[2:n])
hist(dlm$et)
mean(dlm$et[2:n])
ut <- na.omit(dlm$et/sqrt(dlm$Qt))
sd(na.omit(ut))  # it's almost 1


### test data ###

# plot
dlm_test = DLM1U(dECO$p, mu0, C0, bestV, bestW)
sm1U_test = smoother1U(dlm_test)

plot(dECO$p, type="l", ylim = c(0,1))
lines(dlm_test$mt, col = "green", lwd = 2)
lines(sm1U_test[,1], col = "blue", lwd = 2 )
lines((sm1U_test[,1] + 2*sqrt(sm1U_test[,2])), col = "blue", lty=2)
lines((sm1U_test[,1] - 2*sqrt(sm1U_test[,2])), col = "blue", lty=2)

# plot forecast errors
plot(dlm_test$et ~ dECO$t,
     xlab='Year-Month', 
     ylab='Forecast errors',
     main = 'Prevalence data - forecast errors',
     type = "l")
 
# Shewhart control chart
source("EarlyWarning.R")
source("DLM1U.R")
runShewhartUnivariate(dlm_test, 2)

# Plot data and Cusum
cusum = runCusumUnivariate(dlm_test, TRUE)
plot(cusum, type = "l")

# V-mask with reset at alarm
vm2 = runAndPlotV.mask(dlm_test, 8, 0.4, c(-6, 27), TRUE)

##############################################################################
##############################################################################
### DO THE SAME WITH COUNT DATA ###

dlm_c_1 = DLM1U(learn_set$count, mean_count, var_count, 1, 0.04)
sm1U_c_1 = smoother1U(dlm_c_1)

#### improve the function ####

# starting points from smoother function:
mu0_c = sm1U_c_1[1,1]
C0_c = sm1U_c_1[1,2]
    
# Find optimal variance components by EM algorithm
optiEM = runEM1U(learn_set$count, mu0_c, C0_c, 1, 0.04)
bestV_c = optiEM[1]
bestW_c = optiEM[2]

# run filter and smoother with new values
dlm_c_1 = DLM1U(learn_set$count, mu0_c, C0_c, bestV_c, bestW_c)
sm1U_c_1 = smoother1U(dlm_c_1)
plot(learn_set$count, type="l")
lines(dlm_c_1$mt, col = "green", lwd = 2)
lines(sm1U_c_1[,1], col = "blue", lwd = 2 )
lines((sm1U_c_1[,1] + 2*sqrt(sm1U_c_1[,2])), col = "blue", lty=2)
lines((sm1U_c_1[,1] - 2*sqrt(sm1U_c_1[,2])), col = "blue", lty=2)


### assessment ###    

# forecast errors should be normally distributed
n = length(dlm_c_1$et)
plot(dlm_c_1$et[1:(n-1)], dlm_c_1$et[2:n])
hist(dlm_c_1$et)
mean(dlm_c_1$et[2:n])
ut_c <- na.omit(dlm_c_1$et/sqrt(dlm_c_1$Qt))
mean(na.omit(ut_c))
sd(na.omit(ut_c))  # it's almost 1


### test data ###

# plot
dlm_c_test = DLM1U(dECO$count, mu0_c, C0_c, bestV_c, bestW_c)
sm1U_c_test = smoother1U(dlm_c_test)

plot(dECO$count, type="l")
lines(dlm_c_test$mt, col = "green", lwd = 2)
lines(sm1U_c_test[,1], col = "blue", lwd = 2 )
lines((sm1U_c_test[,1] + 2*sqrt(sm1U_c_test[,2])), col = "blue", lty=2)
lines((sm1U_c_test[,1] - 2*sqrt(sm1U_c_test[,2])), col = "blue", lty=2) 

# plot forecast errors
plot(dlm_c_test$et ~ dECO$t,
     xlab='Year-Month', 
     ylab='Forecast errors',
     main = 'Count data - forecast errors',
     type = "l")

# Shewhart control chart
runShewhartUnivariate(dlm_c_test, 2)

# Plot data and Cusum
cusum = runCusumUnivariate(dlm_c_test, TRUE)
plot(cusum, type = "l")

# V-mask with reset at alarm
vm2 = runAndPlotV.mask(dlm_c_test, 8, 0.4, c(-6, 27), TRUE)



