#---------------------------------------------------------------
#---------------------------------------------------------------
#- D A T A     M A N I P U L A T I O N

# This script loads the raw data regarding sugar maple allometries at Hubbard Brook


#- C A L C I U M 
#- foliar Ca content from real field measurements
true.mean = 5
true.std.dev = 0.5
field.sample.size = 12

#- Foliar calcium values are fixed for Monte Carlo, as described in the paper.
#   These are fixed "observations" in the same way that the allometry trees that 
#   Whittaker and colleagues sampled are fixed values.
field.measurements <- rnorm(field.sample.size, true.mean, true.std.dev)
mean(field.measurements)
mean.field <- 5.279
SE.field <- 0.139 
SD.field <- 0.477 





#- D B H    A N D      L E A F      B I O M A S S    
# Whittaker et al. (1974) data on sugar maple leaf and twig biomass as a function of dbh

# Read in the Whittaker et al. (1974) data (sugar maple leaves and twigs)
w.DBH <- c(21.8, 10.7, 15.7, 43.7, 50.8, 28.7, 3.2, 13.9, 23.1, 6.4, 1.9, 29.7, 47, 66)
w.lt <- c(5831, 1130, 1772, 24498, 35468, 10737, 111, 1169, 5401, 310, 78, 8064, 19045, 40002)
wdat <- data.frame(DBH = w.DBH, leaves = w.lt) # combine into dataframe
tree.dbh.for.fake.forest <- sort(w.DBH)[c(2, 4:11, 13)] # sort dbh values in ascending order

# Whittaker et al. (1974) equation parameters for sugar maple leaf and twig biomass 
a = 1.0975
b = 1.9329
whittakersE = 1.385
mn.log.dbh = 1.2449
ssdx.log.dbh = 2.73418
DBH = 17.58

# Fit a linear regression to the log-log data from Whittaker (1974)
summary(lm1 <- lm(log10(w.lt) ~ log10(w.DBH), data = wdat))
wf.a <- as.numeric(coef(lm1)[1]) # parameter 'a' from the model. Note equivalence with hard-coded "a" above.
wf.b <- as.numeric(coef(lm1)[2]) # parameter 'b' from the model. Note equivalence with hard-coded "b" above.
wf.e <- 10^summary(lm1)$sigma[1] # error term from the model
wf.y <- predict(lm1) # fitted values from the model
n <- length(w.lt)    # number of samples
wf.syx <- sqrt( sum(((log10(w.lt)-wf.y)^2)) / (n - 2) ) # sy.x from yanai. Also called the Root Mean Square Error (RMSE). Note log-scale.
denominator <- sum((log10(w.DBH) - mean(log10(w.DBH)))^2) # denominator from eq 5 and 6 in yanai. Note log-scale.

#- values for the analytical solution in Box 3
Xbar <- mean(log10(w.DBH))
ssXa <- sum((log10(w.DBH) - mean(log10(w.DBH)))^2)
sigmax <- sd(log10(w.DBH))

#- end of data manipulation
#---------------------------------------------------------------
#---------------------------------------------------------------
