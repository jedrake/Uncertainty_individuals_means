#---------------------------------------------------------------
#---------------------------------------------------------------
#- M O N T E     C A R L O - (Makes figures 4, 6, and 7)

#- This code utilizes three methods of Monte Carlo uncertainty propagation.
#     Method 1 is based on the uncertainty in the estimation of individuals
#     Method 2 is based on the uncertainty in the estimation of the mean
#     Method 3 is based on the uncertainty in the estimation of the mean, 
#          modified with the uncertainty in individuals. This is the new method
#          that adds a random draw from a distribution determined by the  mean square error
#          of the original allometry to a Monte Carlo based on the mean (i.e., Method 3 = Method 2 + MSE)



# loop through different plot sizes (i.e., different number of observed trees).
# In effect, this loops over the rows in multi-panel figures
for(z in 1: length(plot.sample.size)){ 
  
  no.trees <- plot.sample.size[z] # select the number of individuals ('plot size')
  
  ### Biomass
  dbh <- rep(tree.dbh.for.fake.forest, each = no.trees/10) # selects the same tree dbh values each time. Repeats these values for larger 'plot sizes'.
  
  #- update from Javier on sampling dbh values from a normal distribution
  #dbh <- rep(tree.dbh.for.fake.forest, each = no.trees/10) # selects the same tree dbh values each time. Repeats these values for larger 'plot sizes'.
  #dbh <- rnorm(no.trees,mean(dbh),5)
  
  numerator <- (log10(dbh) - mean(log10(w.DBH)))^2  # How far away is the dbh of each tree from the mean
  #   dbh in Whittaker's dataset?
  
  for(k in 1:no.iter){ # loop over iterations to calculate calcium values for individual trees no.iter times (i.e, 1000 times)
    
    ### Calcium. I am not sure that "calcium3" was calculated properly here. Work on this. 
    calcium1 <- rnorm(plot.sample.size[z], mean.field, SD.field) # for method 1 (individuals)
    calcium2 <- rnorm(1, mean.field, SE.field) # for method 2 (means). Estimate of the population mean for iteration based on the mean and standard error of the field sample 
    #calcium3 <- rnorm(no.trees, calcium2, SD.field) # for method 3 (both). Tree predictions based on the population mean for this iteration and the standard deviation of the field sample
    #calcium3 <- calcium2 + rnorm(no.trees,0,SD.field-SE.field)# for method 3 (both). Tree predictions based on the population mean for this iteration and the standard deviation of the field sample
    calcium3 <- calcium2 + rnorm(no.trees,0,SD.field)# for method 3 (both). Tree predictions based on the population mean for this iteration and the standard deviation of the field sample
    
    
    ### Biomass
    model.error.sample.m2 <- rnorm(1, 0, wf.syx) # for method two mean version
    
    # Randomly sample from the error (individual: sp) for each tree 
    model.error.sample.m1 <- rnorm(no.trees, 0, wf.syx)  # for method one (individuals)
    
    #- calculate error terms as in Yanai et al. 2010 Ecosystems
    sm <- model.error.sample.m2 * sqrt( (1/n) + ( (numerator)/(denominator)  )  ) # yanai eq 5 (mean) [need to log the numerator and demoninator] obtain a different value for each tree
    sp <- model.error.sample.m1 * sqrt( 1 + (1/n) + ( (numerator)/(denominator)  )  ) # yanai eq 6 (individual) [need to log the numerator and demoninator]
    
    # Predict using the Whittaker's equation and the randomly sampled error
    leaf.mass.m0 = 10^(wf.a + wf.b * log10(dbh) )           # without error
    leaf.mass.m1 = 10^(wf.a + wf.b * log10(dbh) + sp)       # using eq 6
    leaf.mass.m2 = 10^(wf.a + wf.b * log10(dbh) + sm)       # using eq 5
    leaf.mass.m3 = 10^(wf.a + wf.b * log10(dbh) + sm + model.error.sample.m1)  # new method #3
    
    # Compute calcium concentrations based on predicted masses and selected calcium
    calcium.11 <- leaf.mass.m1 * calcium1 / 1000000 / 0.002    # convert to kg per ha
    calcium.22 <- leaf.mass.m2 * calcium2 / 1000000 / 0.002
    calcium.33 <- leaf.mass.m3 * calcium3 / 1000000 / 0.002
    
    # Aggregate and output objects
    if(k == 1) { point.c1 <- mean(calcium1)     # Mg/g
    point.c2 <- mean(calcium2)
    point.c3 <- mean(calcium3)
    
    om0.out <- mean(leaf.mass.m0) / 1000 / 0.002 # convert to t per ha
    om1.out <- mean(leaf.mass.m1) / 1000 / 0.002
    om2.out <- mean(leaf.mass.m2) / 1000 / 0.002
    om3.out <- mean(leaf.mass.m3) / 1000 / 0.002
    
    c1.out <- mean(calcium.11)
    c2.out <- mean(calcium.22)
    c3.out <- mean(calcium.33)
    }
    
    if(k > 1 ) { point.c1 <- c(point.c1, mean(calcium1))
    point.c2 <- c(point.c2, mean(calcium2))
    point.c3 <- c(point.c3, mean(calcium3)) 
    
    om0.out <- c(om0.out, mean(leaf.mass.m0) / 1000 / 0.002)
    om1.out <- c(om1.out, mean(leaf.mass.m1) / 1000 / 0.002)
    om2.out <- c(om2.out, mean(leaf.mass.m2) / 1000 / 0.002)
    om3.out <- c(om3.out, mean(leaf.mass.m3) / 1000 / 0.002)
    
    c1.out <- c(c1.out, mean(calcium.11))
    c2.out <- c(c2.out, mean(calcium.22))
    c3.out <- c(c3.out, mean(calcium.33))
    }
    
  } # end k loop over no.iter
  
  #--- put the "data" together
  #--- gdat.z is the focal dataframe that contains all of the estimates of plot-scale sugar maple Ca content
  
  out.z = rbind(data.frame(z = no.trees, method = "Individual", Point.calcium = point.c1, Mass = om1.out, Calcium = c1.out),
                data.frame(z = no.trees, method = "Mean", Point.calcium = point.c2, Mass = om2.out, Calcium = c2.out),
                data.frame(z = no.trees, method = "Both", Point.calcium = point.c3, Mass = om3.out, Calcium = c3.out))
  
  if(z == 1) { gdat.z <- out.z  }
  if(z > 1)  { gdat.z <- rbind(gdat.z, out.z)  }
  
  
} # end z loop plot size

#- define the level of the factor variable for "method"
gdat.z$method <- factor(gdat.z$method,levels=c("Individual","Mean","Both"))


# End Monte Carlo
#---------------------------------------------------------------
#---------------------------------------------------------------











#---------------------------------------------------------------
#---------------------------------------------------------------
# P L O T S    A N D    T A B L E S 

#- Make a series of histograms for leaf mass (This is now Figure 6)
gdat.z %>% 
  mutate(method = factor(method, labels = c("Uncertainty in the prediction\nof an individual", "Uncertainty in the prediction\nof the mean", "Both")), z = factor(z, labels = c("10 trees\n0.002 ha", "30 trees\n0.06 ha", "50 trees\n0.1 ha", "100 trees\n0.2 ha", "1000 trees\n2 ha", "10000 trees\n20 ha"))) %>% 
  ggplot(aes(x = Mass, group = method)) + 
  geom_histogram(binwidth = 200, alpha = 0.4, col = "black") + 
  facet_grid(z ~ method, labeller = label_value) + 
  coord_cartesian(ylim = c(0, 4000), xlim = c(2000, 5800)) + 
  xlab(expression(paste("Leaf mass ", (Mg~ha^-1)))) + 
  ylab("Frequency of Monte Carlo outcomes") + 
  theme_bw()
ggsave("output/Figure6- Regression_graph_leafbiomass_output.png", height = 8, width = 7)
#ggsave("output/Figure6- Regression_graph_leafbiomass_output.tiff", units="in", width=8, height=7, dpi=300)


#- Make a series of histograms for leaf Ca concentration. This is now Figure 4.
gdat.z %>% 
  mutate(method = factor(method, labels = c("Uncertainty in the prediction\nof an individual", "Uncertainty in the prediction\nof the mean", "Both")), z = factor(z, labels = c("10 trees\n0.002 ha", "30 trees\n0.06 ha", "50 trees\n0.1 ha", "100 trees\n0.2 ha", "1000 trees\n2 ha", "10000 trees\n20 ha"))) %>% 
  ggplot(aes(x = Point.calcium, group = method)) + 
  geom_histogram(binwidth = 0.1, alpha = 0.4, col = "black") + 
  facet_grid(z ~ method, labeller = label_value) + 
  coord_cartesian(ylim = c(0, 3000), xlim = c(4.5, 6)) + 
  xlab(expression(paste("Leaf calcium concentration ", (mg~g^-1)))) + 
  ylab("Frequency of Monte Carlo outcomes") + 
  theme_bw()
ggsave("output/Figure4- Regression_graph_calcium_concentration_output.png", height = 8, width = 7)
#ggsave("output/Figure4- Regression_graph_calcium_concentration_output.tiff", units="in", width=8, height=7, dpi=300)

#- Make a series of histograms for total leaf Ca content. This is now Figure 7.
gdat.z %>% 
  mutate(method = factor(method, labels = c("Uncertainty in the prediction\nof an individual", "Uncertainty in the prediction\nof the mean", "Both")), z = factor(z, labels = c("10 trees\n0.002 ha", "30 trees\n0.06 ha", "50 trees\n0.1 ha", "100 trees\n0.2 ha", "1000 trees\n2 ha", "10000 trees\n20 ha"))) %>% 
  ggplot(aes(x = Calcium, group = method)) + 
  geom_histogram(binwidth = 2, alpha = 0.4, col = "black") + 
  facet_grid(z ~ method, labeller = label_value) + 
  coord_cartesian(ylim = c(0, 5000), xlim = c(5, 35)) + 
  xlab(expression(paste("Leaf calcium content ", (kg~ha^-1)))) + 
  ylab("Frequency of Monte Carlo outcomes") + 
  theme_bw()
ggsave("output/Figure7- Regression_graph_calcium_content_output.png", height = 8, width = 7)
#ggsave("output/Figure7- Regression_graph_calcium_content_output.tiff", units="in", width=8, height=7, dpi=300)

# Write out a Table of SD values from Monte Carlo output values
gdat.z %>% 
  group_by(method, z) %>% 
  summarise_all(list(sd=sd)) %>% 
  write.table("output/SD_values_Monte_Carlo.csv", sep = ",", row.names = F)
#---------------------------------------------------------------
#---------------------------------------------------------------
