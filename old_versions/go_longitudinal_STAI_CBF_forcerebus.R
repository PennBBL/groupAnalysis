##############################################################################
################                                               ###############
################   GO1-GO2 Voxelwise Analysis for CBF          ###############
################           Angel Garcia de la Garza            ###############
################              angelgar@upenn.edu               ###############
################                 03/01/2016                    ###############
##############################################################################

##############################################################################
################                  Load Libraries               ###############
##############################################################################

library(ggplot2)
library(base)
library(reshape2)
library(gamm4)
library(stats)
library(knitr)
library(mgcv)
library(visreg)
library(plyr)
library(ANTsR)
library(parallel)


##############################################################################
################        Define variable paths                  ###############
##############################################################################

subjDataName<-"/home/angelgar/STAI_CBF_Voxelwise/input/final_voxelwise_cbf_stai_longitudinal.Rds"  #subject data frame 
InDirRoot<-"/home/angelgar/STAI_CBF_Voxelwise/input"  #input directory
OutDirRoot<-"/home/angelgar/STAI_CBF_Voxelwise/output"  #output directory
outName<-"testrun_forted" #name of output subfolder

##############################################################################
################        Find path to files                    ###############
##############################################################################
#assumes naming from flameo
imageName<-paste(InDirRoot,"allgolongitudinalcbf_4mm.nii.gz",sep="/")  #this assumes that you have concatenated all of the images into a 4d-timeseries.   
# **WARNING** must be careful and have order here the same as teh subject list and data or else will permute your data.  can make this with fslmerge -t
# would suggest downsampling this 4D timeseries to 4mm NOT 2mm for speed--- can use fslmaths [inputTimesereis] -subsamp2 [output timeseries]
imageIds<-paste(InDirRoot,"go_longitudinal_final_paths.csv",sep="/")  #this is list of paths to scans must be in same order as our 4D Image
maskName<-paste(InDirRoot,"n404_longitudinal_cbf_gm_mask_final2.nii.gz",sep="/")  #this is a mask to run analyses within.  I will copy a file for you


##############################################################################
################            Load subject data                  ###############
##############################################################################
subjData<-readRDS(subjDataName) ##Read Data
subjs<-read.csv(imageIds, header=F, as.is=T) ##This is a list of paths, you should also have this same list on the data


#Assumes data is in the same order as paths you will double check this 
#check versus subject list
for(i in 1:dim(subjData)[1]) {
  if(!identical(subjs$V1[i],subjData$path[i])){
    print(i)
    print("subjectIds not matched-- ERROR!!!")
    stop()
  }
}




##############################################################################
################        Make Output Directory                  ###############
##############################################################################
outDir<-paste(OutDirRoot,outName,sep="/")
logDir<-paste(OutDirRoot,"logs",sep="/")
#Will return a warning if logDir and outDir already exist
dir.create(outDir)
dir.create(logDir)

##############################################################################
################              Echo Arguments                   ###############
##############################################################################
print("Arguments are:")
print(paste("subject data is:",subjDataName))
print(paste("input images are:",imageName))
print(paste("mask is:",maskName))
print(paste("ouput directory is:",outDir))
print(paste("log directory is:",logDir))


###cleanup logdir
system(paste('rm -f', file.path(logDir, '*')))


##############################################################################
################                   Load data                   ###############
##############################################################################
mask<-antsImageRead(maskName,3)
imageIn<-antsImageRead(imageName,4)
imageMat<-timeseries2matrix(imageIn,mask)


##############################################################################
################           Preallocate output                  ###############
##############################################################################
pOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)
tOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)


##############################################################################
################        Run model using parallel               ###############
##############################################################################
timeOn<-proc.time()
covariates=" ~ stai_stai_tr +sex +s(age) +s(age, by=sex)"


# We create a list of formulas for each voxel in our data. 
# Each element in the list will have formula with a different voxel as the dependent variable
m <- mclapply(1:dim(imageMat)[2], function(x) {as.formula(paste(paste0("imageMat[,",x,"]"), covariates, sep=""))})

# We run gamm4 for each element in our list of formulas
# This will return a list with each element being a table of p-values
# This can be customized based on what you need form the model 
# Never save the whole gamm4 object it will be too big and will mclapply
model <- mclapply(m, function(x) {summary(gamm4(formula = x, data=subjData, REML=T, random = ~ (1|bblid))$gam)$p.table}, mc.cores = 10)

loopTime<-proc.time()-timeOn


##############################################################################
################        Allocate out t-map and z-map            ###############
##############################################################################


for(i in 1:length(model)){
  pOut[i,1]<-model[[i]][which(rownames(model[[i]]) == "stai_stai_tr"),4]
  tOut[i,1]<-model[[i]][which(rownames(model[[i]]) == "stai_stai_tr"),3]
}



##############################################################################
################             Write out maps                    ###############
##############################################################################
pOutImage<-antsImageClone(mask)
tOutImage<-antsImageClone(mask)
pOutImage[mask==1]<-pOut
tOutImage[mask==1]<-tOut

setwd(outDir)
antsImageWrite(pOutImage,"gamm4P.nii.gz")
antsImageWrite(tOutImage,"gamm4T.nii.gz")
