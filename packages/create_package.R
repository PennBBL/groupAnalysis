install.packages("devtools")
install.packages("roxygen2")

setwd("~/scripts/rPackages/")
create("voxel")

setwd("~/scripts/rPackages/voxel/R")


library(devtools)
library(roxygen2)

setwd("~/scripts/rPackages/voxel/")
document()

setwd("~/scripts/rPackages/")
install("voxel")
library(voxel)

library(mgcv)
library(oro.nifti)
library(parallel)
library(gamm4)

image <- nifti(img = array(1:1600, dim =c(4,4,4,25)))
mask <- nifti(img = array(1:4, dim = c(4,4,4,1)))
set.seed(1)
covs <- data.frame(x = runif(25), id = rep(1:5,5))
formula <- "~ s(x)"
randomFormula <- "~(1|id)"
models <- gamCluster(image=image, mask=mask, covsFormula=formula, subjData=covs, ncores = 1, REML=T)

