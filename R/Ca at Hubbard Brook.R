#---------------------------------------------------------------
#---------------------------------------------------------------
#  M A K E    F I N A L   F I G U R E   -  H U B B A R D      B R O O K

#- Take the excel sheet that was developed by many people in Ruth's lab
#   and most recently by the student Joe Nash. This excel sheet is included as
#   a supplementary material to the manuscript and called XXXX. 
#   Do the Monte Carlo simulation there, extract the data (column AD of any of the three
#    "iterations" tabs, which has total Ca in Kg ha-1), and graph the
#    distributions. I pasted the total ecosystem Ca data 
#    into a tab labeled "total Ca data" and saved that tab as an csv file.

joedat <- read.csv("data/joenashdata_10000.csv")

#- define the level of the factor variable for "method"
joedat$Method <- factor(joedat$Method,levels=c("Individual","Mean","Both"))


#- Make a series of histograms for total Ca content

#- set up pdf for the plot
pdf(file="output/Figure7- Ca at HB.pdf",width=8,height=6)
layout(matrix(1:3, 1, 3, byrow = TRUE), 
       widths=c(1,1,1), heights=c(1,1,1))
par(mar=c(0,2,1,0),oma=c(7,9,5,5),cex.lab=1.7,cex.axis=1.2)

xlims=c(550,800)
ylims=c(0,3500)
binwidth = 20

hist(subset(joedat,Method=="Individual")$Ca_kg_ha,main="",
     xlim=xlims,ylim=ylims,freq=T,
     breaks = seq(from=530, to=850, by=5))
hist(subset(joedat,Method=="Mean")$Ca_kg_ha,main="",
     xlim=xlims,ylim=ylims,
     breaks = seq(from=530, to=850, by=binwidth))
hist(subset(joedat,Method=="Both")$Ca_kg_ha,main="",
     xlim=xlims,ylim=ylims,
     breaks = seq(from=530, to=850, by=binwidth))

#- add top labels
mtext("Individuals",side=3,line=0,cex=1.3,adj=0.1,outer=T)
mtext("Mean",side=3,line=0,cex=1.3,adj=0.5,outer=T)
mtext("Both",side=3,line=0,cex=1.3,adj=0.92,outer=T)

#- add axis labels
mtext(expression(Total~Ca~content~"("*kg~ha^-1*")"),side=1,line=4,cex=1.3,outer=T)
mtext("Frequency of Monte Carlo outcomes (10000 simulations)",side=2,line=3,cex=1.3,outer=T)


dev.off()
#---------------------------------------------------------------
#---------------------------------------------------------------




#---------------------------------------------------------------
#---------------------------------------------------------------
#- calculate the CV's based on Joe's MonteCarlo analysis for the text
cv.dat <- summaryBy(Ca_kg_ha~Method,data=joedat,FUN=c(mean,sd))
cv.dat$CV <- with(cv.dat,Ca_kg_ha.sd/Ca_kg_ha.mean)
cv.dat
#---------------------------------------------------------------
#---------------------------------------------------------------
