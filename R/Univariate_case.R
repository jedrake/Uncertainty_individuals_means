#---------------------------------------------------------------
#---------------------------------------------------------------
# This script assesses the univariate case of the Ca concentration of leaves
#  in sugar maple trees. It creates Figure 2 for the manuscript.
#---------------------------------------------------------------
#---------------------------------------------------------------



#---------------------------------------------------------------
#---------------------------------------------------------------
#  M A K E    F I G U R E    2 - C A L C I U M     V E R S I O N

#- Illustrate propagating the uncertainty of means, individuals, and both.
#

# define a distribution based on a real set of measurements. N=12 trees
mean.field <- 5.279
SE.field <- 0.139 
SD.field <- 0.477

trees_vec <- c(5,10,30,50,100,1000,10000) # number of individuals in each iteration

#- in the uncertainty of individuals, we will randomly sample from a normal distribution defined by
#   a mean of 5.279 and the SD.field of 0.477 FOR EACH TREE SEPARATELY

#- in the uncertainty of means, we will randomly from a normal distribution defined by
#   a mean of 5.279 and the SE.field of 0.139. his would be done once per each iteration and applied to ALL trees. 

#  In the "both" case, simulate the uncertainty of means by randomly drawing from a normal distribution defined by
#   a mean of 5.279 and the SE.field of 0.139. This will be a random estimate of the mean based on the sample, that 
#   will be applied to all individual trees in that iteration. Then, there is a second random draw from a 
#   normal distribution defined by this random estimate of the mean based on the sample (for that iteration), 
#   with a variance term defined as the difference between 
#   SD.field (0.477) and SE.field (0.139), or 0.338. This could be done by adding rnorm(mean=0,se=SD.field-SE.field)

# loop over "trees"
out_calcium1 <- out_calcium2 <- out_calcium3 <- list()
for (a in 1:length(trees_vec)){
  
  trees <- trees_vec[a]
  # loop over iterations
  calcium1 <- calcium2 <- calcium3 <- list()
  for (i in 1:no.iter){
    
    calcium1[[i]] <- rnorm(trees, mean.field, SD.field) # for method 1 (individuals)
    calcium2[[i]] <- rep(rnorm(1, mean.field, SE.field),trees) # for method 2 (means). Estimate of the population mean for iteration based on the mean and standard error of the field sample 
    #calcium3[[i]] <- calcium2[[i]] + rnorm(trees, 0, SD.field^2-SE.field^2) # for method 3 (both). Tree predictions based on the population mean for this iteration and the standard deviation of the field sample
    calcium3[[i]] <- calcium2[[i]] + rnorm(trees, 0, SD.field) # for method 3 (both). Tree predictions based on the population mean for this iteration and the standard deviation of the field sample
  }
  #- extract dataframe from list, calculate mean of each iteration
  calcium1.df <- as.data.frame(do.call(rbind,calcium1)) 
  calcium1.df$mean <- rowMeans(calcium1.df[,1:trees])
  calcium1.df$trees <- trees
  
  calcium2.df <- as.data.frame(do.call(rbind,calcium2))
  calcium2.df$mean <- rowMeans(calcium2.df[,1:trees])
  calcium2.df$trees <- trees
  
  calcium3.df <- as.data.frame(do.call(rbind,calcium3))
  calcium3.df$mean <- rowMeans(calcium3.df[,1:trees])
  calcium3.df$trees <- trees
  
  out_calcium1[[a]] <- calcium1.df[,c("trees","mean")]
  out_calcium2[[a]] <- calcium2.df[,c("trees","mean")]
  out_calcium3[[a]] <- calcium3.df[,c("trees","mean")]
}
calcium1.df <- as.data.frame(do.call(rbind,out_calcium1)) 
calcium2.df <- as.data.frame(do.call(rbind,out_calcium2)) 
calcium3.df <- as.data.frame(do.call(rbind,out_calcium3)) 

variances1 <- summaryBy(mean~trees,data=calcium1.df,FUN=var)
variances2 <- summaryBy(mean~trees,data=calcium2.df,FUN=var)
variances3 <- summaryBy(mean~trees,data=calcium3.df,FUN=var)

#- what should we have gotten?
variances2$mean.var*(1-1/variances2$trees)

#- So we varied the number of trees in our new sample (k) to be 10, 30, 50, 100, 1000, and 10000.
#  And we calculated the variances of the three approaches
variances1
variances2
variances3
0.477^2*(1/12+1/c(5, 10, 30, 50, 100, 1000, 10000))


#- set up pdf for the plot. 
pdf(file="output/Figure2- Univariate case.pdf",width=8,height=8)
layout(matrix(1:9, 3, 3, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1,1))
par(mar=c(0,0,0,0),oma=c(7,9,5,5),cex.lab=1.7,cex.axis=1.2)


ylims <- c(4.3,6.3)
symbol <- 16

#- plot the three methods for the first set of tree numbers (5 trees)
plot(subset(calcium1.df,trees==trees_vec[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
abline(h=mean.field,col="darkgrey")

plot(subset(calcium2.df,trees==trees_vec[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean.field,col="darkgrey")

plot(subset(calcium3.df,trees==trees_vec[1])$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean.field,col="darkgrey")

#- plot the three methods for 30 trees
plot(subset(calcium1.df,trees==30)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
abline(h=mean.field,col="darkgrey")

plot(subset(calcium2.df,trees==30)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean.field,col="darkgrey")

plot(subset(calcium3.df,trees==30)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=F,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean.field,col="darkgrey")

#- plot the three methods for 10000 trees
plot(subset(calcium1.df,trees==10000)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=T,tcl=0.5)
abline(h=mean.field,col="darkgrey")

plot(subset(calcium2.df,trees==10000)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean.field,col="darkgrey")

plot(subset(calcium3.df,trees==10000)$mean,ylim=ylims,xaxt="n",yaxt="n", frame.plot=TRUE, xlab="", ylab="",pch=symbol,col=alpha("black",0.25))
axis(side=1,labels=T,tcl=0.5);axis(side=2,labels=F,tcl=0.5)
abline(h=mean.field,col="darkgrey")


#- add titles
mtext("Individuals",side=3,line=0,cex=1.5,adj=0.07,outer=T)
mtext("Mean",side=3,line=0,cex=1.5,adj=0.5,outer=T)
mtext("Combined",side=3,line=0,cex=1.5,adj=0.92,outer=T)

mtext("10 Trees",side=2,line=3,cex=1,adj=0.9,outer=T)
mtext("30 Trees",side=2,line=3,cex=1,adj=0.5,outer=T)
mtext("10000 Trees",side=2,line=3,cex=1,adj=0.1,outer=T)

mtext("Iterations",side=1,line=3,cex=1.5,adj=0.5,outer=T)

mtext(expression(Ca~concentration~(mg~g^-1)),side=2,line=5,cex=1.5,adj=0.5,outer=T)


dev.off()

# E N D    O F    U N I V A R I A T E     C A S E    (FIGURE 2)
#---------------------------------------------------------------
#---------------------------------------------------------------
