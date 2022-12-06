

#---------------------------------------------------------------
#---------------------------------------------------------------
#  M A K E    F I G U R E    1 - color option

# I decided to drop this version, as the code would rescale the color ramp for each panel, which I didn't like. 

#- set the color palette
# greys.ramp <- colorRampPalette(c("whitesmoke","blue", "black"),alpha=0.7)
# 
# 
# 
# #- set up pdf for the plot. Work on making this more colorful/transparent to address overplotting
# pdf(file="output/Figure1_color.pdf",width=8,height=8)
# layout(matrix(1:9, 3, 3, byrow = TRUE), 
#        widths=c(1,1,1), heights=c(1,1,1))
# par(mar=c(0,0,0,0),oma=c(7,9,5,5),cex.lab=1.7,cex.axis=1.2)
# 
# 
# ylims <- c(4.3,6.3)
# symbol <- 16
# 
# #- plot the three methods for the first set of tree numbers (5 trees)
# 
# #- create color densities
# dc1 <- densCols(1:no.iter,subset(calcium1.df,trees==trees_vec[1])$mean,
#                   colramp=greys.ramp)
# dc2 <- densCols(1:no.iter,subset(calcium2.df,trees==trees_vec[1])$mean,
#                 colramp=greys.ramp)
# dc3 <- densCols(1:no.iter,subset(calcium3.df,trees==trees_vec[1])$mean,
#                 colramp=greys.ramp)
# dc4 <- densCols(1:no.iter,subset(calcium1.df,trees==30)$mean,
#                 colramp=greys.ramp)
# dc5 <- densCols(1:no.iter,subset(calcium2.df,trees==30)$mean,
#                 colramp=greys.ramp)
# dc6 <- densCols(1:no.iter,subset(calcium3.df,trees==30)$mean,
#                colramp=greys.ramp)
# dc7 <- densCols(1:no.iter,subset(calcium1.df,trees==10000)$mean,
#                 colramp=greys.ramp)
# dc8 <- densCols(1:no.iter,subset(calcium2.df,trees==10000)$mean,
#                 colramp=greys.ramp)
# dc9 <- densCols(1:no.iter,subset(calcium3.df,trees==10000)$mean,
#                 colramp=greys.ramp)
# 
# plot(subset(calcium1.df,trees==trees_vec[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc1)
# axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# plot(subset(calcium2.df,trees==trees_vec[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc2)
# axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# plot(subset(calcium3.df,trees==trees_vec[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc3)
# axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# #- plot the three methods for 30 trees
# plot(subset(calcium1.df,trees==30)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc4)
# axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# plot(subset(calcium2.df,trees==30)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc5)
# axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# plot(subset(calcium3.df,trees==30)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc6)
# axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# #- plot the three methods for 10000 trees
# plot(subset(calcium1.df,trees==10000)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc7)
# axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# plot(subset(calcium2.df,trees==10000)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc8)
# axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# plot(subset(calcium3.df,trees==10000)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=dc8)
# axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
# abline(h=mean.field,col="darkgrey")
# 
# 
# #- add titles
# mtext("Individuals",side=3,line=0,cex=1.5,adj=0.07,outer=T)
# mtext("Mean",side=3,line=0,cex=1.5,adj=0.5,outer=T)
# mtext("Combined",side=3,line=0,cex=1.5,adj=0.92,outer=T)
# 
# mtext("5 Trees",side=2,line=3,cex=1,adj=0.9,outer=T)
# mtext("30 Trees",side=2,line=3,cex=1,adj=0.5,outer=T)
# mtext("10000 Trees",side=2,line=3,cex=1,adj=0.1,outer=T)
# 
# mtext("Iterations",side=1,line=3,cex=1.5,adj=0.5,outer=T)
# 
# mtext(expression(Ca~concentration~(mg~g^-1)),side=2,line=5,cex=1.5,adj=0.5,outer=T)
# 
# 
# dev.off()

#---------------------------------------------------------------
#---------------------------------------------------------------








#---------------------------------------------------------------
#---------------------------------------------------------------
#  M A K E    F I G U R E    1 - M A S S   V E R S I O N 

# This has been removed from the manuscript.

#- Illustrate propagating the uncertainty of means, individuals, and both.
#

# define a distribution based on a real set of measurements
wdat
wdat$DBH_log <- log10(wdat$DBH)
wdat$leaves_log <- log10(wdat$leaves)


#- get randomly samples values of the slope and intercept
lm.leaves <- lm(leaves_log~DBH_log,data=wdat)
a <-    1.09752
b <-    1.93294
Sigma <- vcov(lm.leaves)
coeffs_set <- mvrnorm(n=100, c(a, b), Sigma)             # random sampling pairs of a & b from multiple normal distribution


no.trees <- c(10,30,10000) # number of individuals in each iteration

no.iter = 1000

# loop over "trees"
out_mass1 <- out_mass2 <- out_mass3 <- list()
for (t in 1:length(no.trees)){
  
  trees <- no.trees[t]
  # loop over iterations
  dbh <- rep(tree.dbh.for.fake.forest, each = trees/10) # selects the same tree dbh values each time. Repeats these values for larger 'plot sizes'.
  
  mass1 <- mass2 <- mass3 <- list()
  for (i in 1:no.iter){
    
    #- get separate estimates of the slope and intercept for each tree
    coeffs_ind <- mvrnorm(n=length(dbh), c(a, b), Sigma)            
    
    #- get a single estimate of the slope and intercept for this iteration
    coeffs_mean <- mvrnorm(n=1, c(a, b), Sigma)            
    
    
    #- calculate leaf biomass
    mass1[[i]] = 10^(coeffs_ind[,1] + coeffs_ind[,2] * log10(dbh) + rnorm(length(dbh),0,wf.syx)) # for individuals, note the additional uncertainty of the RMSE
    mass2[[i]] = 10^(coeffs_mean[1] + coeffs_mean[2] * log10(dbh) )           
    mass3[[i]] = 10^(coeffs_mean[1] + coeffs_mean[2] * log10(dbh) + rnorm(length(dbh),0,wf.syx))           
    
  }
  #- extract dataframe from list, calculate mean of each iteration
  mass1.df <- as.data.frame(do.call(rbind,mass1)) 
  mass1.df$mean <- rowMeans(mass1.df[,1:trees])
  mass1.df$trees <- trees
  
  mass2.df <- as.data.frame(do.call(rbind,mass2))
  mass2.df$mean <- rowMeans(mass2.df[,1:trees])
  mass2.df$trees <- trees
  
  mass3.df <- as.data.frame(do.call(rbind,mass3))
  mass3.df$mean <- rowMeans(mass3.df[,1:trees])
  mass3.df$trees <- trees
  
  out_mass1[[t]] <- mass1.df[,c("trees","mean")]
  out_mass2[[t]] <- mass2.df[,c("trees","mean")]
  out_mass3[[t]] <- mass3.df[,c("trees","mean")]
}
mass1.df <- as.data.frame(do.call(rbind,out_mass1)) 
mass2.df <- as.data.frame(do.call(rbind,out_mass2)) 
mass3.df <- as.data.frame(do.call(rbind,out_mass3)) 




#- set up pdf for the plot
pdf(file="output/Figure1-mass.pdf",width=8,height=8)
layout(matrix(1:9, 3, 3, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1,1))
par(mar=c(0,0,0,0),oma=c(7,7,5,5),cex.lab=1.7,cex.axis=1.2)


ylims <- c(3000,13000)
symbol <- 3

#- plot the three methods for the first set of tree numbers
plot(subset(mass1.df,trees==no.trees[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
abline(h=mean(subset(mass1.df,trees==no.trees[1])$mean),col="darkgrey")

plot(subset(mass2.df,trees==no.trees[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean(subset(mass2.df,trees==no.trees[1])$mean),col="darkgrey")

plot(subset(mass3.df,trees==no.trees[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean(subset(mass3.df,trees==no.trees[1])$mean),col="darkgrey")

#- plot the three "mean" methods
plot(subset(mass1.df,trees==no.trees[2])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
abline(h=mean(subset(mass1.df,trees==no.trees[2])$mean),col="darkgrey")

plot(subset(mass2.df,trees==no.trees[2])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean(subset(mass2.df,trees==no.trees[2])$mean),col="darkgrey")

plot(subset(mass3.df,trees==no.trees[2])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean(subset(mass3.df,trees==no.trees[2])$mean),col="darkgrey")

#- plot the three "both" methods
plot(subset(mass1.df,trees==no.trees[3])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
abline(h=mean(subset(mass3.df,trees==no.trees[3])$mean),col="darkgrey")

plot(subset(mass2.df,trees==no.trees[3])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean(subset(mass3.df,trees==no.trees[3])$mean),col="darkgrey")

plot(subset(mass3.df,trees==no.trees[3])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol)
axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean(subset(mass3.df,trees==no.trees[3])$mean),col="darkgrey")


#- add titles
mtext("Individuals",side=3,line=0,cex=1.5,adj=0.07,outer=T)
mtext("Mean",side=3,line=0,cex=1.5,adj=0.5,outer=T)
mtext("Combined",side=3,line=0,cex=1.5,adj=0.92,outer=T)

mtext("10 Trees",side=2,line=3,cex=1,adj=0.9,outer=T)
mtext("30 Trees",side=2,line=3,cex=1,adj=0.5,outer=T)
mtext("10000 Trees",side=2,line=3,cex=1,adj=0.1,outer=T)

mtext("Iterations",side=1,line=3,cex=1.5,adj=0.5,outer=T)


dev.off()

# E N D    O F    F I G U R E    1- mass version (removed from manuscript)
#---------------------------------------------------------------
#---------------------------------------------------------------






# 
# 
# #---------------------------------------------------------------
# #---------------------------------------------------------------
# # demonstrate that 95% of observations fall within prediction interval
# 
# #- generate a fake forest of trees
# dbh <- runif(10000,min=0.2,max=1.5) # uniform distribution from 0.2 to 1.5 (dbh on log scale)
# leaves <- 1.09752+1.93294*dbh + rnorm(10000,mean=0,sd=0.3)
# 
# fakelm <- lm(leaves~dbh)
# plot(leaves~dbh)
# 
# z <- predict(fakelm, newdata=data.frame(dbh), se.fit=TRUE)
# 
# se.CI <- z$se.fit
# se.PI <- sqrt(z$se.fit^2 + z$residual.scale^2) # not sure where this comes from. See https://stackoverflow.com/questions/38109501/how-does-predict-lm-compute-confidence-interval-and-prediction-interval
# 
# 
# # compute confidence/prediction intervals at 95% level
# alpha <- qt((1-0.95)/2, df = z$df) # note that df is n-2
# CI_upr <- z$fit + -alpha * se.CI
# CI_lwr <- z$fit + alpha * se.CI
# PI_upr <- z$fit + -alpha * se.PI
# PI_lwr <- z$fit + alpha * se.PI
# 
# fakedat <- data.frame(dbh=dbh,leaves=leaves,CI_upr,CI_lwr,PI_upr,PI_lwr)
# fakedat$within_CI <- ifelse(fakedat$leaves>fakedat$CI_lwr & fakedat$leaves<fakedat$CI_upr,1,0)
# fakedat$within_PI <- ifelse(fakedat$leaves>fakedat$PI_lwr & fakedat$leaves<fakedat$PI_upr,1,0)
# 
# # 94.99% of observations are within the prediction interval
# table(fakedat$within_PI)
# 0.95*10000
# 
# # 0.2% of observations are within the confidence interval
# table(fakedat$within_CI)
# 
# #- pull out data that are within 0.05 of 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, and 1.4 (dbh log units)
# subset(fakedat,round(dbh))
#---------------------------------------------------------------
#---------------------------------------------------------------






























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

#- end of data manipulation
#---------------------------------------------------------------
#---------------------------------------------------------------


















# 
# #---------------------------------------------------------------
# #---------------------------------------------------------------
# #- M O N T E     C A R L O ,  U S I N G     S L O P E    A N D    I N T E R C E P T
# 
# #- This code utilizes three methods of monte carlo uncertainty propogation.
# #     Method 1 applies random slope and intercept terms for each tree (prediction of individual)
# #     Method 2 applies random slope and intercept terms for each iteration (prediction of mean)
# #     Method 3 is ___
# 
# 
# 
# # loop through different plot sizes (i.e., different number of observed trees).
# # In effect, this loops over the rows in multi-panel figures
# for(z in 1: length(plot.sample.size)){ 
#   
#   no.trees <- plot.sample.size[z] # select the number of individuals ('plot size')
#   
#   ### Biomass
#   dbh <- rep(tree.dbh.for.fake.forest, each = no.trees/10) # selects the same tree dbh values each time. Repeats these values for larger 'plot sizes'.
#   
# 
#   
#   for(k in 1:no.iter){ # loop over interations to calculate calcium values for individual trees no.iter times (i.e, 1000 times)
#     
#     
#     
#     #- get separate estimates of the slope and intercept for each tree
#     coeffs_ind <- mvrnorm(n=length(dbh), c(a, b), Sigma)            
#     
#     #- get a single estimate of the slope and intercept for this iteration
#     coeffs_mean <- mvrnorm(n=1, c(a, b), Sigma)            
#     
#     
#     #- calculate leaf biomass
#     leaf.mass.m0 = 10^(coeffs_ind[,1] + coeffs_ind[,2] * log10(dbh) + rnorm(length(dbh),0,wf.syx)) # for individuals, note the additional uncertainty of the RMSE
#     leaf.mass.m1 = 10^(coeffs_mean[1] + coeffs_mean[2] * log10(dbh) )           
#     leaf.mass.m2 = 10^(coeffs_mean[1] + coeffs_mean[2] * log10(dbh) + rnorm(length(dbh),0,wf.syx))           
#     
#     # 
#     # ### Calcium. I am not sure that "calcium3" was calculated properly here. Work on this. 
#     # calcium1 <- rnorm(plot.sample.size[z], mean.field, SD.field) # for method 1 (individuals)
#     # calcium2 <- rnorm(1, mean.field, SE.field) # for method 2 (means). Estimate of the population mean for iteration based on the mean and standard error of the field sample 
#     # #calcium3 <- rnorm(no.trees, calcium2, SD.field) # for method 3 (both). Tree predictions based on the population mean for this iteration and the standard deviation of the field sample
#     # calcium3 <- calcium2 + rnorm(no.trees,0,SD.field-SE.field)# for method 3 (both). Tree predictions based on the population mean for this iteration and the standard deviation of the field sample
#     # 
#     # 
#     # 
#     # ### Biomass
#     # model.error.sample.m2 <- rnorm(1, 0, wf.syx) # for method two mean version
#     # 
#     # # Randomly sample from the error (individual: sp) for each tree 
#     # model.error.sample.m1 <- rnorm(no.trees, 0, wf.syx)  # for method one (individuals)
#     # 
#     # #- calculate error terms as in Yanai et al. 2010 Ecosystems
#     # sm <- model.error.sample.m2 * sqrt( (1/n) + ( (numerator)/(denominator)  )  ) # yanai eq 5 (mean) [need to log the numerator and demoninator] obtain a different value for each tree
#     # sp <- model.error.sample.m1 * sqrt( 1 + (1/n) + ( (numerator)/(denominator)  )  ) # yanai eq 6 (individual) [need to log the numerator and demoninator]
#     # 
#     # # Predict using the Whittaker's equation and the randomly sampled error
#     # leaf.mass.m0 = 10^(wf.a + wf.b * log10(dbh) )           # without error
#     # leaf.mass.m1 = 10^(wf.a + wf.b * log10(dbh) + sp)       # using eq 6
#     # leaf.mass.m2 = 10^(wf.a + wf.b * log10(dbh) + sm)       # using eq 5
#     # leaf.mass.m3 = 10^(wf.a + wf.b * log10(dbh) + sm + model.error.sample.m1)  # new method #3
#     # 
#     # # Compute calcium concentrations based on predicted masses and selected calcium
#     # calcium.11 <- leaf.mass.m1 * calcium1 / 1000000 / 0.002    # convert to kg per ha
#     # calcium.22 <- leaf.mass.m2 * calcium2 / 1000000 / 0.002
#     # calcium.33 <- leaf.mass.m3 * calcium3 / 1000000 / 0.002
#     # 
#     # Aggregate and output objects
#     if(k == 1) { 
#     #point.c1 <- mean(calcium1)     # Mg/g
#     #point.c2 <- mean(calcium2)
#     #point.c3 <- mean(calcium3)
#     
#     om0.out <- mean(leaf.mass.m0) / 1000 / 0.002 # convert to t per ha
#     om1.out <- mean(leaf.mass.m1) / 1000 / 0.002
#     om2.out <- mean(leaf.mass.m2) / 1000 / 0.002
#     #om3.out <- mean(leaf.mass.m3) / 1000 / 0.002
#     
#     #c1.out <- mean(calcium.11)
#     #c2.out <- mean(calcium.22)
#     #c3.out <- mean(calcium.33)
#     }
#     
#     if(k > 1 ) { 
#     #point.c1 <- c(
#     #point.c1, mean(calcium1))
#     #point.c2 <- c(point.c2, mean(calcium2))
#     #point.c3 <- c(point.c3, mean(calcium3)) 
#     
#     om0.out <- c(om0.out, mean(leaf.mass.m0) / 1000 / 0.002)
#     om1.out <- c(om1.out, mean(leaf.mass.m1) / 1000 / 0.002)
#     om2.out <- c(om2.out, mean(leaf.mass.m2) / 1000 / 0.002)
#     #om3.out <- c(om3.out, mean(leaf.mass.m3) / 1000 / 0.002)
#     
#     #c1.out <- c(c1.out, mean(calcium.11))
#     #c2.out <- c(c2.out, mean(calcium.22))
#     #c3.out <- c(c3.out, mean(calcium.33))
#     }
#     
#   } # end k loop over no.iter
#   
#   #--- put the "data" together
#   #--- gdat.z is the focal dataframe that contains all of the estimates of plot-scale sugar maple Ca content
#   
#   
#   out.z = rbind(data.frame(z = no.trees, method = "Individual", Mass = om0.out),
#                               data.frame(z = no.trees, method = "Mean", Mass = om1.out),
#                               data.frame(z = no.trees, method = "Both", Mass = om2.out))
#                 
#   #out.z = rbind(data.frame(z = no.trees, method = "Individual", Point.calcium = point.c1, Mass = om1.out, Calcium = c1.out),
#   #              data.frame(z = no.trees, method = "Mean", Point.calcium = point.c2, Mass = om2.out, Calcium = c2.out),
#   #              data.frame(z = no.trees, method = "Both", Point.calcium = point.c3, Mass = om3.out, Calcium = c3.out))
#   
#   if(z == 1) { gdat.z <- out.z  }
#   if(z > 1)  { gdat.z <- rbind(gdat.z, out.z)  }
#   
#   
# } # end z loop plot size
# 
# 
# # End Monte Carlo
# #---------------------------------------------------------------
# #---------------------------------------------------------------
# 
# 
# #---------------------------------------------------------------
# #---------------------------------------------------------------
# # P L O T S    A N D    T A B L E S 
# 
# #- Make a series of histograms for leaf mass
# gdat.z %>% 
#   mutate(method = factor(method, labels = c("Uncertainty in the prediction\nof an individual", "Uncertainty in the prediction\nof the mean","Both")), z = factor(z, labels = c("10 trees\n0.002 ha", "30 trees\n0.06 ha", "50 trees\n0.1 ha", "100 trees\n0.2 ha", "1000 trees\n2 ha", "10000 trees\n20 ha"))) %>% 
#   ggplot(aes(x = Mass, group = method)) + 
#   geom_histogram(binwidth = 200, alpha = 0.4, col = "black") + 
#   facet_grid(z ~ method, labeller = label_value) + 
#   coord_cartesian(ylim = c(0, 4000), xlim = c(2000, 6000)) + 
#   xlab(expression(paste("Leaf mass ", (Mg~ha^-1)))) + 
#   ylab("Frequency of Monte Carlo outcomes") + 
#   theme_bw()
# ggsave("Regression_graph_biomass_output.png", height = 8, width = 7)
