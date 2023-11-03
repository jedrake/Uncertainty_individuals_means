#---------------------------------------------------------------
#---------------------------------------------------------------
#- This script illustrates the allometry and makes figure 4
#---------------------------------------------------------------
#---------------------------------------------------------------



#---------------------------------------------------------------
#---------------------------------------------------------------
#  M A K E    F I G U R E    4 - ( illustrates allometry )




#   This is meant to be a demonstration of this approach graphically,
#     in the context of regression confidence and prediction intervals.
#   It uses the Whittaker data that were read in previously by allometric_data.R
wdat
wdat$DBH_log <- log10(wdat$DBH)
wdat$leaves_log <- log10(wdat$leaves)



#- get randomly samples values of the slope and intercept
lm.leaves <- lm(leaves_log~DBH_log,data=wdat)
a <-    1.09752
b <-    1.93294
Sigma <- vcov(lm.leaves)
coeffs_set <- mvrnorm(n=100, c(a, b), Sigma)             # random sampling pairs of a & b from multiple normal distribution


# calculate error in the prediciton of the mean and individual
denominator <- sum((wdat$DBH_log-mean(wdat$DBH_log))^2)
nrow(wdat) # 14 observations
MSE <- sum((wdat$leaves_log-wf.y)^2) / (14-2) # should the denominator be n-1 or n-2
VarYhat_mean <- MSE * (1/14 + ((wdat$DBH_log-mean(wdat$DBH_log))^2 / denominator  ))
VarYhat_ind <- MSE * (1 + 1/14 + ((wdat$DBH_log-mean(wdat$DBH_log))^2 / denominator  ))

z <- predict(lm1, newdata=data.frame(wdat$DBH_log), se.fit=TRUE)


# calculate the confidence interval (from VarYhat_mean) and the prediction interval (from VarYhat_ind)
alpha <- qt((1-0.95)/2, df = z$df) # note that df is n-2
ymax_mean <- wf.y + -alpha * sqrt(VarYhat_mean)
ymin_mean <- wf.y + alpha * sqrt(VarYhat_mean)
ymax_ind <- wf.y + -alpha * sqrt(VarYhat_ind)
ymin_ind <- wf.y + alpha * sqrt(VarYhat_ind)




#--------------------------------------------------------
# extract confidence and prediction interval via normal R tools

#- the "by hand" method below exactly recomputes the base R tools.
CI_r <- data.frame(predict(lm1,newdata=data.frame(wdat$DBH_log),interval="confidence"))
PI_r <- data.frame(predict(lm1,newdata=data.frame(wdat$DBH_log),interval="prediction"))


se.CI <- z$se.fit
se.PI <- sqrt(z$se.fit^2 + z$residual.scale^2) # not sure where this comes from. See https://stackoverflow.com/questions/38109501/how-does-predict-lm-compute-confidence-interval-and-prediction-interval


# compute confidence/prediction intervals at 95% level
alpha <- qt((1-0.95)/2, df = z$df) # note that df is n-2
CI_upr <- z$fit + -alpha * se.CI
CI_lwr <- z$fit + alpha * se.CI
PI_upr <- z$fit + -alpha * se.PI
PI_lwr <- z$fit + alpha * se.PI

#--------------------------------------------------------


#- Now make a graph with three areas


#- function to add polygons
addpoly <- function(x,y1,y2,col=alpha("lightgrey",0.7),...){
  ii <- order(x)
  y1 <- y1[ii]
  y2 <- y2[ii]
  x <- x[ii]
  polygon(c(x,rev(x)), c(y1, rev(y2)), col=col, border=NA,...)
}



#- illustrate the uncertainty in mean and individual predictions
#windows(100,100)
pdf(file="output/Figure5- allometry.pdf",width=8,height=8)
#tiff(file="output/Figure5- allometry.tiff", units="in", width=8, height=8, res=300)


par(mar=c(6,6,1,1),cex.axis=2)
plot(leaves_log~DBH_log,data=wdat,pch=1,ylim=c(1,5),xlab="",ylab="")

#- add lines for slope and intercept
#for (i in 1:nrow(coeffs_set)){
#  abline(a=coeffs_set[i,1],b=coeffs_set[i,2],col="grey")
#}
points(leaves_log~DBH_log,data=wdat,pch=15,ylim=c(1,5))

title(xlab=expression(DBH~(log[10]~(cm))),cex.lab=2.5,line=4)
title(ylab=expression(Leaf~biomass~(log[10]~(g))),cex.lab=2.5,line=3)

#- add polygons for uncertainty in individual (red) and mean (blue)
addpoly(x=wdat$DBH_log, y1=ymax_ind,y2=ymin_ind, col=alpha("lightgoldenrod1",0.25))
addpoly(x=wdat$DBH_log, y1=ymax_mean,y2=ymin_mean, col=alpha("lightsalmon",0.25))
lines(wf.y~wdat$DBH_log,lty=1,ylim=c(1,5),xlab="",ylab="",cex=1.5)
points(leaves_log~DBH_log,data=wdat,pch=16,cex=2)
legend("topleft",legend=c("Prediction","Confidence"),fill=c("lightgoldenrod1","lightsalmon"),cex=2)


#- re-scale to linear, make plot
ymax_ind_linear <- 10^ymax_ind
ymin_ind_linear <- 10^ymin_ind
ymax_mean_linear <- 10^ymax_mean
ymin_mean_linear <- 10^ymin_mean
pred_linear <- 10^wf.y
#---------------------------------------------------------------
#---------------------------------------------------------------

dev.off()
#- E N D     OF     F I G U R E     5  (demonstration of allometry)
#---------------------------------------------------------------
#---------------------------------------------------------------
