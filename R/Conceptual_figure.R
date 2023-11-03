#---------------------------------------------------------------
#---------------------------------------------------------------
#- This script makes a conceptual figure showing a distribution
#   and the samples that we took from that distribution. 
#   It makes Figure 1.
#---------------------------------------------------------------
#---------------------------------------------------------------





#---------------------------------------------------------------
#---------------------------------------------------------------
#  M A K E    C O N C E P T U A L   F I G U R E -  D I S T R I B U T I O N S



#------ Plot the "true" mean
#Create a sequence of 100 equally spaced numbers between -4 and 4
x <- seq(2, 8, length=100)

#create a vector of values that shows the height of the probability distribution
#for each value in x
y <- dnorm(x,mean=5,sd=0.5)




#plot x and y as a scatterplot with connected lines (type = "l") and add
#an x-axis with custom labels
pdf(file="output/Figure1- distributions.pdf",width=8,height=8)
#tiff(file="output/Figure1- distributions.tiff", units="in", width=8, height=8, res=300)

plot(x,y, type = "l", lwd = 2, axes = FALSE, 
     xlab = (expression(paste("Leaf calcium concentration ", (mg~g^-1)))), ylab = "")
axis(1, at = 0:8, labels = seq(0,8))
text(x=5,y=0.85, expression(paste("True mean, ", mu, "=5")),xpd=T)
abline(v=5,lty=2)


#----- overlay an illustration of the sampled samples
samples <- rnorm(12,mean=5.279, sd =  0.477)
points(samples,rep(0.01,length(samples)))
est.mean <- mean(samples)
est.sd <- sd(samples)
est.se <- sd(samples/sqrt(length(samples)))
points(est.mean,0.01,pch=16,cex=1.5)
#adderrorbars(x=est.mean,y=0.03,SE=est.se*1.96,direction="leftright")
text(x=5.1,y=0.048, expression(paste("Estimated mean, ", bar("x"), " =5.3")),xpd=T)


dev.off()
#---------------------------------------------------------------
#---------------------------------------------------------------
