

#- create output directory, if it doesn't exist
if(!dir.exists("Output")) dir.create(file.path("Output"),showWarnings=F)



#- load libraries
Library <- function(pkg, ...){
  
  PACK <- .packages(all.available=TRUE)
  pkgc <- deparse(substitute(pkg))
  
  if(pkgc %in% PACK){
    library(pkgc, character.only=TRUE)
  } else {
    install.packages(pkgc, ...)
    library(pkgc, character.only=TRUE)
  }
  
}

#- load some libraries
Library(ggplot2)
Library(dplyr)
Library(scales)
Library(MASS)
Library(doBy)
Library(grDevices)
Library(doBy)
source("R/generic_functions.R")
