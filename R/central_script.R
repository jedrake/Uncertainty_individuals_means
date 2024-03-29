#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------
#- This script demonstrates three ways of propagating the uncertainty of total foliar Ca of 
#     forested stands that contain different numbers of individual trees.
#
#  This is the central script for the by Yanai et al. (2023) in Ecosystems, entitled:
#   "Propagating uncertainty in predicting individuals and means illustrated with foliar chemistry and forest biomass"
#
#  This script calls many other scripts that perform subroutines, as described below. By default,
#   the scripts called here will generate pdfs for figures in the "output/" folder. If you prefer,
#   you can have the scripts output high-resolution TIFFs (300 dpi) by uncommenting the tiff() line and commenting
#   the pdf() line in each script. Note that the 300 dpi TIFFs are large files (>10 MB each). 
#  
#  This code was primarily written by Hannah Buckley in 2016 and was edited by John Drake in 2018-2023.
#---------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------



#---------------------------------------------------------------
#---------------------------------------------------------------
#- load some libraries that are required
source("R/load_libraries.R")
#---------------------------------------------------------------
#---------------------------------------------------------------



#---------------------------------------------------------------
#---------------------------------------------------------------
#- P A R A M E T E R S 
#- Set the number of iterations for Monte Carlo loop
no.iter = 10000

#- Set the number of trees in each simulations. Simulates small to very large forest plots
plot.sample.size = c(10, 30, 50, 100, 1000, 10000) # units are trees per plot

#- end of parameters
#---------------------------------------------------------------
#---------------------------------------------------------------



#---------------------------------------------------------------
#---------------------------------------------------------------
#- Allometric data. Load and manipulate the underlying raw
#   allometric data for sugar maple trees at Hubbard Brook
source("R/allometric_data.R")
#---------------------------------------------------------------
#---------------------------------------------------------------


#---------------------------------------------------------------
#---------------------------------------------------------------
#- Make the conceptual figure showing the distribution. Figure 1.
source("R/Conceptual_figure.R")
#---------------------------------------------------------------
#---------------------------------------------------------------



#---------------------------------------------------------------
#---------------------------------------------------------------
#- Assess the univariate case and make Figure 3. Note that Figure 2
#   is a flow-chart diagram that is not represented in this code.

#  For illustration purposes, we use fewer iterations here than
#   in subsequent analyses for which we advocate using 10,000 iterations.
no.iter <- 3000
source("R/Univariate_case.R")
no.iter <- 10000 # set the number of iterations back to 10000 for subsequent scripts
#---------------------------------------------------------------
#---------------------------------------------------------------



#---------------------------------------------------------------
#---------------------------------------------------------------
#- Illustrate the allometry and make Figure 5.
source("R/illustrate_allometry.R")
#---------------------------------------------------------------
#---------------------------------------------------------------



#---------------------------------------------------------------
#---------------------------------------------------------------
#- Run the main Monte Carlo uncertainty assessment code.
#  This one may take some time to run (~10 min) in a way that 
#  depends strongly on the number of iterations (no.iter).

#  Create Figures 4, 6, and 7.
source("R/main_monte_carlo.R") 
#---------------------------------------------------------------
#---------------------------------------------------------------




#---------------------------------------------------------------
#---------------------------------------------------------------
#- Run an alternative version of the main Monte Carlo code,
#  but in an assessment of the analytical solution presented in 
#  box 3

source("R/main_monte_carlo_analytical_assessment.R") 
#---------------------------------------------------------------
#---------------------------------------------------------------


#---------------------------------------------------------------
#---------------------------------------------------------------
#- Make the last figure regarding total Ca content at Hubbard Brook.
#  The code to do the simulation is embedded in an excel sheet.
#  This script plots the results. Creates Figure 8
source("R/Ca at Hubbard Brook.R")
#---------------------------------------------------------------
#---------------------------------------------------------------