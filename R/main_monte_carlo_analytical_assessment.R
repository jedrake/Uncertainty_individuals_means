#---------------------------------------------------------------
#---------------------------------------------------------------
#- M O N T E     C A R L O - (Makes figures 3, 5, and 6)

#- This code utilizes three methods of Monte Carlo uncertainty propagation.
#     Method 1 is based on the uncertainty in the estimation of individuals
#     Method 2 is based on the uncertainty in the estimation of the mean
#     Method 3 is based on the uncertainty in the estimation of the mean, 
#          modified with the uncertainty in individuals. This is the new method
#          that adds a random draw from a distribution determined by the  mean square error
#          of the original allometry to a Monte Carlo based on the mean (i.e., Method 3 = Method 2 + MSE)


library(dplyr)
# loop through different plot sizes (i.e., different number of observed trees).
# In effect, this loops over the rows in multi-panel figures
for(z in 1: length(plot.sample.size)){ 
  #k
  no.trees <- plot.sample.size[z] # select the number of individuals ('plot size')
  
  ### Biomass
  ##ATTENTION!! Play between line 24 (original, used throughout) 
  #             and 25 (pure sampling, no repetitions, used to assess analytical solution in box 3)
  #dbh <- rep(tree.dbh.for.fake.forest, each = no.trees/10) # selects the same tree dbh values each time. Repeats these values for larger 'plot sizes'.
  dbh <- sample(tree.dbh.for.fake.forest, 10,TRUE)
  dbh <- rnorm(no.trees,mean(dbh),2)
  
  numerator <- (log10(dbh) - mean(log10(w.DBH)))^2  # How far away is the dbh of each tree from the mean

  #For assessment of analytical solution of box 3
  var_x=var(log10(dbh))##This is SIGMA_X (in log(dbh))
  meanlogdbh<-mean(log10(dbh))
  var=wf.syx^2 #model variance
  n=n #no. of sample trees
  k=no.trees
  SSz=denominator
  beta1=wf.b
  meanz=mean(log10(w.DBH))
  mu=meanlogdbh
  #End assessment of analytical solution of box 3
  
  for(iter in 1:no.iter){ # loop over iterations to calculate calcium values for individual trees no.iter times (i.e, 1000 times)
    ### Calcium.  
    calcium1 <- rnorm(plot.sample.size[z], mean.field, SD.field) # for method 1 (individuals)
    calcium2 <- rnorm(1, mean.field, SE.field) # for method 2 (means). Estimate of the population mean for iteration based on the mean and standard error of the field sample 
    calcium3 <- calcium2 + rnorm(no.trees,0,SD.field)# for method 3 (both). Tree predictions based on the population mean for this iteration and the standard deviation of the field sample
    
    ### Biomass
    # Randomly sample from the error (individual: sp) for each tree 
    model.error.sample.m1 <- rnorm(no.trees, 0, wf.syx)  # for method one (individuals), wf.syx is the sigma of the mdoel
    mean(model.error.sample.m1)
    
    # Randomly sample from the error (mean: sm) for each tree 
    model.error.sample.m2 <- rnorm(1, 0, wf.syx) # for method two mean version, wf.syx is the sigma of the mdoel

    #- calculate error terms as in Yanai et al. 2010 Ecosystems
    sm <- model.error.sample.m2 * sqrt( (1/n) + ( (numerator)/(denominator)  )  ) # yanai eq 5 (mean) [need to log the numerator and demoninator] obtain a different value for each tree
    sp <- model.error.sample.m1 * sqrt( 1 + (1/n) + ( (numerator)/(denominator)  )  ) # yanai eq 6 (individual) [need to log the numerator and demoninator]
    
    # Predict using the Whittaker's equation and the randomly sampled error
    #leaf.mass.m0 = 10^(wf.a + wf.b * log10(dbh) )           # without error
    leaf.mass.m1 = 10^(wf.a + wf.b * log10(dbh) + sp)       # using eq 6: individuals
    leaf.mass.m2 = 10^(wf.a + wf.b * log10(dbh) + sm)       # using eq 5: mean
    leaf.mass.m3 = 10^(wf.a + wf.b * log10(dbh) + sm + model.error.sample.m1)  #method #3: 
    
    #Specifically on a log-scale, without back transformation
    #Logleaf.mass.m0 = wf.a + wf.b * log10(dbh)            # without error
    Logleaf.mass.m1 = wf.a + wf.b * log10(dbh) + sp       # using eq 6: individuals
    Logleaf.mass.m2 = wf.a + wf.b * log10(dbh) + sm      # using eq 5: mean
    Logleaf.mass.m3 = wf.a + wf.b * log10(dbh) + sm + model.error.sample.m1  # new method #3
    
    
    # Compute calcium concentrations based on predicted masses and selected calcium
    #calcium.00 <- leaf.mass.m0 * calcium0 / 1000000 / 0.002    # convert to kg per ha. Individuals
    calcium.11 <- leaf.mass.m1 * calcium1 / 1000000 / 0.002    # convert to kg per ha. Individuals
    calcium.22 <- leaf.mass.m2 * calcium2 / 1000000 / 0.002   # convert to kg per ha. Mean
    calcium.33 <- leaf.mass.m3 * calcium3 / 1000000 / 0.002   # convert to kg per ha. BOTH
    
    # Aggregate and output objects
    if(iter == 1) {
      point.c1 <- mean(calcium1)     # Mg/g. Individuals
    point.c2 <- mean(calcium2)                    # Mg/g.  Mean
    point.c3 <- mean(calcium3)                    # Mg/g.  BOTH
    
    lm1.out <- mean(leaf.mass.m1) # kg, Individuals
    lm2.out <- mean(leaf.mass.m2)  # kg, Mean
    lm3.out <- mean(leaf.mass.m3)  # kg, BOTH
    
    Loglm1.out <- mean(Logleaf.mass.m1) # kg, Individuals
    Loglm2.out <- mean(Logleaf.mass.m2)  # kg, Mean
    Loglm3.out <- mean(Logleaf.mass.m3)  # kg, BOTH
    
    
    om1.out <- mean(leaf.mass.m1) / 1000 / 0.002 # convert to t per ha, Individuals
    om2.out <- mean(leaf.mass.m2) / 1000 / 0.002 # convert to t per ha, Mean
    om3.out <- mean(leaf.mass.m3) / 1000 / 0.002 # convert to t per ha, BOTH
    
    c1.out <- mean(calcium.11)
    c2.out <- mean(calcium.22)
    c3.out <- mean(calcium.33)
    }
    
    if(iter > 1 ) {
    point.c1 <- c(point.c1, mean(calcium1)) #Individuals
    point.c2 <- c(point.c2, mean(calcium2))               #Mean
    point.c3 <- c(point.c3, mean(calcium3))               #BOTH
    
    
    lm1.out <- c(lm1.out, mean(leaf.mass.m1)) #Biomass, kg. Individuals
    lm2.out <- c(lm2.out, mean(leaf.mass.m2)) #Biomass, kg. Mean
    lm3.out <- c(lm3.out, mean(leaf.mass.m3)) #Biomass, kg. BOTH
    
    Loglm1.out <- c(Loglm1.out, mean(Logleaf.mass.m1)) #Biomass, kg. Individuals
    Loglm2.out <- c(Loglm2.out, mean(Logleaf.mass.m2)) #Biomass, kg. Mean
    Loglm3.out <- c(Loglm3.out, mean(Logleaf.mass.m3)) #Biomass, kg. BOTH
    
    
    om1.out <- c(om1.out, mean(leaf.mass.m1) / 1000 / 0.002) #Biomass, T/ha. Individuals
    om2.out <- c(om2.out, mean(leaf.mass.m2) / 1000 / 0.002) #Biomass, T/ha. Mean
    om3.out <- c(om3.out, mean(leaf.mass.m3) / 1000 / 0.002) #Biomass, T/ha. BOTH
    
    c1.out <- c(c1.out, mean(calcium.11))   #Calcium Kg/ha . Individuals
    c2.out <- c(c2.out, mean(calcium.22))   #Calcium Kg/ha . Mean
    c3.out <- c(c3.out, mean(calcium.33))   #Calcium Kg/ha . BOTH
    }
    
  } # end iter loop over no.iter
  
  #--- put the "data" together
  #--- gdat.z is the focal dataframe that contains all of the estimates of plot-scale sugar maple Ca content
  
  out.z = rbind(data.frame(k = no.trees, method = "Individual", Point.calcium = point.c1, leafmass=lm1.out, Logleafmass=Loglm1.out, Mass = om1.out, Calcium = c1.out),
                data.frame(k = no.trees, method = "Mean", Point.calcium = point.c2, leafmass=lm2.out,Logleafmass=Loglm2.out,Mass = om2.out, Calcium = c2.out),
                data.frame(k = no.trees, method = "Both", Point.calcium = point.c3, leafmass=lm3.out,Logleafmass=Loglm3.out,Mass = om3.out, Calcium = c3.out))
  
  if(z == 1) { gdat.z <- out.z  }
  if(z > 1)  { gdat.z <- rbind(gdat.z, out.z)  }
  
  #Analytical result: Terry's development in box 3
  ana_result.z=var/n+var*var_x/(no.trees*SSz)+var_x*beta1^2/no.trees+var*(mu-meanz)^2/SSz+var/no.trees
  
  outana.z=rbind(data.frame(k=no.trees,
                             first=var/n, 
                             second=var*(mu-meanz)^2/SSz,
                             third=var*var_x/(no.trees*SSz),
                             fourth=var_x*beta1^2/no.trees,
                             fifth=var/no.trees,
                             total=ana_result.z))
  
  if(z == 1) { gdatana.z <- outana.z  }
  if(z > 1)  { gdatana.z <- rbind(gdatana.z, outana.z)  }  
  
    
    
} # end z loop plot size

#- define the level of the factor variable for "method"
gdat.z$method <- factor(gdat.z$method,levels=c("Individual","Mean","Both"))

# End Monte Carlo
#---------------------------------------------------------------
#---------------------------------------------------------------



#Here we compared both regression results, the MC (gdat.z) and the analytical (gdatanal.z).
#We compare Log (leaf mass)
mcresults <- gdat.z %>% group_by(method,k) %>% summarise(var=var(Logleafmass, na.rm = TRUE))
gdatana.z

plot(subset(mcresults,method=="Both")$var~gdatana.z$total,
     xlab="Analytical solution (variance)",
     ylab="Monte Carlo solution (variance)")
abline(0,1)
