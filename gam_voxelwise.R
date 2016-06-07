##############################################################################
################                                               ###############
################   GO1-GO2 Voxelwise Analysis for CBF          ###############
################           Angel Garcia de la Garza            ###############
################              angelgar@upenn.edu               ###############
################                 05/10/2016                    ###############
##############################################################################

suppressMessages(require(optparse))

##############################################################################
################                 Option List                   ###############
##############################################################################


option_list = list(
  make_option(c("-c", "--covariates"), action="store", default=NA, type='character',
              help="Full path to RDS covariate file.  
              Please include bblid and scanid in this file as a column each"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="Full path to output directory"), 
  make_option(c("-p", "--imagepaths"), action="store", default=NA, type='character',
              help="Name of the variable in the covariate file that contains the path to the images to be analyzed"), 
  make_option(c("-m", "--mask"), action="store", default=NA, type='character',
              help="Full path to mask"), 
  make_option(c("-s", "--smoothing"), action="store", default=NA, type='numeric',
              help="The smoothing in sigmas required for the fourd image. Please write in 0 if no smoothing is wanted"), 
  make_option(c("-i", "--inclusion"), action="store", default=NA, type='character',
              help="Name of inclusion variable on dataset. By default 1 means include. This will subset your rds file"), 
  make_option(c("-u", "--subjId"), action="store", default=NA, type='character',
              help="subjID name on the covariates dataset"), 
  make_option(c("-f", "--formula"), action="store", default=NA, type='character',
              help="Formula for covariates to be used, should only include the right hand side of the formula.
              DO NOT INCLUDE SPACES IN THE FORMULA
              Example: ~ stai_stai_tr+sex+s(age)+s(age,by=sex)"), 
  make_option(c("-n", "--numbercores"), action="store", default=10, type='numeric',
              help="Number of cores to be used, default is 10")
  
  )

opt = parse_args(OptionParser(option_list=option_list))

for (i in 1:length(opt)){
  if (is.na(opt)[i] == T) {
    cat('User did not specify all arguments.\n')
    cat('Use gam_voxelwise.R -h for an expanded usage menu.\n')
    quit()
  }
}


##############################################################################
################                  Load Libraries               ###############
##############################################################################

suppressMessages(require(ggplot2))
suppressMessages(require(base))
suppressMessages(require(reshape2))
suppressMessages(require(nlme))
suppressMessages(require(lme4))
suppressMessages(require(gamm4))
suppressMessages(require(stats))
suppressMessages(require(knitr))
suppressMessages(require(mgcv))
suppressMessages(require(plyr))
suppressMessages(require(ANTsR))
suppressMessages(require(parallel))
suppressMessages(require(optparse))


##############################################################################
################              Declare Variables               ###############
##############################################################################

subjDataName <- opt$covariates
OutDirRoot <- opt$output
namePaths <- opt$imagepaths
maskName  <- opt$mask
smooth <- opt$smoothing
inclusionName <- opt$inclusion
subjID <- opt$subjId
covsFormula <- opt$formula
randomFormula <- opt$random
ncores <- opt$numbercores

print("All arguments have been read")

##############################################################################
################            Load subject data                  ###############
##############################################################################
subjData<-readRDS(subjDataName) ##Read Data
subset <- which(subjData[inclusionName] == 1) ##Find subset for analysis
subjData <- subjData[subset, ] #subset data

print("Covariates file has been read")

##############################################################################
################    Create Analysis Directory                  ###############
##############################################################################

OutDir <- paste0(OutDirRoot, "/n",dim(subjData)[1],"_",namePaths,"_",inclusionName,"_smooth",as.character(smooth))
dir.create(OutDir)
setwd(OutDir)

##############################################################################
################     Create and output fourd image             ###############
##############################################################################

fourdcommand <- subjData[1, namePaths]
for (i in 2:dim(subjData)[1]) {
  fourdcommand <- paste(fourdcommand, subjData[i, namePaths])
}

fourdcommand <- paste('fslmerge -t fourd.nii.gz', fourdcommand)

system(fourdcommand, wait=T)
print("4d image succesfully created")

system(paste0("scp ", maskName," ",OutDir, "/mask.nii.gz"), wait=T)
print("mask succesfully copied")

if (smooth > 0) {
  system(paste0("fslmaths ",OutDir, "/fourd.nii.gz -s ",smooth," ",OutDir,"/fourd.nii.gz"), wait=T)
  print("Smoothing of 4d image done")
} else {
  print("No smoothing done")
}




##############################################################################
################            Output summary files               ###############
##############################################################################

write.csv(subjData[, namePaths], paste0(namePaths,".csv"), row.names = F)
write.csv(subjData[, subjID], paste0(subjID,".csv"), row.names = F)

print("Succesfully wrote paths and id files")

##############################################################################
################        Make Output Directory                  ###############
##############################################################################

outName <- gsub("~", "", covsFormula)
outName <- gsub("\\+","_",outName)
outName <- gsub("\\(","",outName)
outName <- gsub("\\)","",outName)
outName <- gsub(",","",outName)
outName <- gsub("=","",outName)
outName <- gsub("\\*","and",outName)
outName <- gsub(":","and",outName)

random <- gsub("~", "", randomFormula)
random <- gsub("\\(", "", random)
random <- gsub("\\)", "", random)
random <- gsub("\\|", "", random)

outsubDir <- paste0("n",dim(subjData)[1],"gam_Cov_",outName)

outsubDir<-paste(OutDir,outsubDir,sep="/")

logDir<-paste(OutDir,"logs",sep="/")

#Will return a warning if logDir and outsubDir already exist
dir.create(logDir)
dir.create(outsubDir)

system(paste('rm -f', file.path(outsubDir, '*')))

##############################################################################
################              Echo Arguments                   ###############
##############################################################################
system( paste0("echo Arguments are: >> ", outsubDir,"/logs.txt"))
system( paste0("echo Covariates file is: ", subjDataName,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Output directory is: ", OutDir,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Path name variable in covarites file is: ", namePaths,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Mask path is: ", maskName,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Smoothing is: ", smooth," >> ", outsubDir,"/logs.txt"))
system( paste0("echo Inclusion variable name is: ", inclusionName,">> ", outsubDir,"/logs.txt"))
system( paste0("echo ID variable name is: ", subjID,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Formula for fixed effects is: ", outName,">> ", outsubDir,"/logs.txt"))
system( paste0("echo Number of cores is: ", ncores," >> ", outsubDir,"/logs.txt"))


###cleanup logdir
system(paste('rm -f', file.path(logDir, '*')))

##############################################################################
################                   Load data                   ###############
##############################################################################

maskName <- paste0(OutDir,"/mask.nii.gz")
imageName <- paste0(OutDir,"/fourd.nii.gz")


mask<-antsImageRead(maskName,3)
imageIn<-antsImageRead(imageName,4)
imageMat<-timeseries2matrix(imageIn,mask)

print("Data has been loaded")

##############################################################################
################           Preallocate output                  ###############
##############################################################################
pOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)
tOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)
zOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)

print("Preallocate output done")

##############################################################################
################        Run model using parallel               ###############
##############################################################################
timeOn<-proc.time()

length.voxel <- ceiling(dim(imageMat)[2]/20)

# We create a list of formulas for each voxel in our data. 
# Each element in the list will have formula with a different voxel as the dependent variable

m <- mclapply(1:10, function(x) {as.formula(paste(paste0("imageMat[,",x,"]"), covsFormula, sep=""))}, mc.cores = ncores)
model <- mclapply(m, function(x) {
  foo <- summary(gam(formula = x, data=subjData, method="REML"))
  return(rbind(foo$p.table,foo$s.table))
}, mc.cores = ncores)


for (k in 1:20) {
  if (k == 20) {
  m <- mclapply((11 + (k-1)*length.voxel):dim(imageMat)[2], function(x) {as.formula(paste(paste0("imageMat[,",x,"]"), covsFormula, sep=""))}, mc.cores = ncores)  
  } else {
  m <- mclapply((11 + (k-1)*length.voxel):(10 + (k)*length.voxel), function(x) {as.formula(paste(paste0("imageMat[,",x,"]"), covsFormula, sep=""))}, mc.cores = ncores)  
  }
  model.temp <- mclapply(m, function(x) {
    foo <- summary(gam(formula = x, data=subjData, method="REML"))
    return(rbind(foo$p.table,foo$s.table))
  }, mc.cores = ncores)
  
  model <- c(model, model.temp)
}

loopTime<-proc.time()-timeOn


print("Models are done")
print(loopTime/60)



##############################################################################
################        Allocate out t-map and z-map            ###############
##############################################################################

setwd(outsubDir)

for (j in 1:dim(model[[1]])[1]) {
  
  variable <- rownames(model[[1]])[j]
  
  if (grepl("s(", variable, fixed=T)) {
    
    for(i in 1:length(model)){
      pOut[i,1]<- model[[i]][which(rownames(model[[i]]) == variable),4]
      zOut[i,1]<- qnorm(model[[i]][which(rownames(model[[i]]) == variable),4], lower.tail=F)
    }
    
    pOutImage<-antsImageClone(mask)
    pOutImage[mask==1]<-pOut
    
    zOutImage<-antsImageClone(mask)
    zOutImage[mask==1]<-zOut
    
    var <- gsub("\\(", "", variable)
    var <- gsub("\\)", "", var)
    var <- gsub(",", "", var)
    var <- gsub("=", "", var)
    
    
    
    antsImageWrite(pOutImage,paste0("gamP_",var,".nii.gz"))
    antsImageWrite(zOutImage,paste0("gamZ_",var,".nii.gz"))
  }
  else {
    for(i in 1:length(model)){
      pOut[i,1]<-model[[i]][which(rownames(model[[i]]) == variable),4]
      zOut[i,1]<-qnorm(model[[i]][which(rownames(model[[i]]) == variable),4], lower.tail=F)
      tOut[i,1]<-model[[i]][which(rownames(model[[i]]) == variable),3]
    }
    
    pOutImage<-antsImageClone(mask)
    zOutImage<-antsImageClone(mask)
    tOutImage<-antsImageClone(mask)
    pOutImage[mask==1]<-pOut
    zOutImage[mask==1]<-zOut
    tOutImage[mask==1]<-tOut
    
    var <- gsub("\\(", "", variable)
    var <- gsub("\\)", "", var)
    var <- gsub("\\*","and",var)
    var <- gsub(":","and",var)
    
    antsImageWrite(pOutImage,paste0("gamP_",var,".nii.gz"))
    antsImageWrite(zOutImage,paste0("gamZ_",var,".nii.gz"))
    antsImageWrite(tOutImage,paste0("gamT_",var,".nii.gz"))
    
  }
}

print("Write t-maps and p-maps is done")
