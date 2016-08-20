##############################################################################
################                                               ###############
################             GAM Voxelwise Wrapper             ###############
################           Angel Garcia de la Garza            ###############
################              angelgar@upenn.edu               ###############
################                 07/11/2016                    ###############
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
              Example: ~ stai_stai_tr+sex+s(age)+s(age,by=sex)"),
  make_option(c("-a", "--padjust"), action="store", default="none", type='character',
              help="method used to adjust pvalues, default is `none`"),
  make_option(c("-k", "--splits"), action="store", default=10, type='numeric',
              help="number of splits to divide the data in, default is 10. To minimize data usage"),
  make_option(c("-n", "--numbercores"), action="store", default=10, type='numeric',
              help="Number of cores to be used, default is 10"),
  make_option(c("-d", "--skipfourD"), action="store", default=FALSE, type='logical',
              help="Option to skip creation of fourdD image and look for it in the Analysis Directory.
              4D image must be labeled as 'fourd.nii.gz'. Will also skip smoothing step.
              Default (FALSE) means to not skip")
  )

opt = parse_args(OptionParser(option_list=option_list))

for (i in 1:length(opt)){
  if (is.na(opt)[i] == T) {
    cat('User did not specify all arguments.\n')
    cat('Use gam_voxelwise.R -h for an expanded usage menu.\n')
    quit()
  }
}

print("##############################################################################")
print("################  Generalized Additive Model Voxelwise Script  ###############")
print("################            Angel Garcia de la Garza           ###############")
print("################              angelgar@upenn.edu               ###############")
print("################                 Version 2.0.0                 ###############")
print("##############################################################################")

##############################################################################
################                  Load Libraries               ###############
##############################################################################

print("Loading Libraries")

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
suppressMessages(require(oro.nifti))
suppressMessages(require(parallel))
suppressMessages(require(optparse))
suppressMessages(require(fslr))



##############################################################################
################              Declare Variables               ###############
##############################################################################

print("Reading Arguments")

subjDataName <- opt$covariates
OutDirRoot <- opt$output
namePaths <- opt$imagepaths
maskName  <- opt$mask
smooth <- opt$smoothing
inclusionName <- opt$inclusion
subjID <- opt$subjId
covsFormula <- opt$formula
pAdjustMethod <- opt$padjust
splits <- opt$splits
ncores <- opt$numbercores
skipFourD <- opt$skipfourD

methods <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY","fdr", "none")
if (!any(pAdjustMethod == methods)) {
  print("p.adjust.method is not a valid one, reverting back to 'none'")
  pAdjustMethod <- "none"
}


##############################################################################
################            Load subject data                  ###############
##############################################################################

print("Loading covariates file")
subjData<-readRDS(subjDataName) ##Read Data
subset <- which(subjData[inclusionName] == 1) ##Find subset for analysis
subjData <- subjData[subset, ] #subset data



##############################################################################
################    Create Analysis Directory                  ###############
##############################################################################


print("Creating Analysis Directory")
OutDir <- paste0(OutDirRoot, "/n",dim(subjData)[1],"_",namePaths,"_",inclusionName,"_smooth",as.character(smooth))
dir.create(OutDir)
setwd(OutDir)

##############################################################################
################     Create and output fourd image             ###############
##############################################################################


if (!skipFourD) {
  print("Merging images and saving out a fourd image")
  subjList <- as.character(subjData[,grep(namePaths, names(subjData))])
  length.subj <- length(subjList)
  k <- splits
  break.subj <- ceiling(length.subj / k)
  
  subMergeNames <- "foo"
  for (i in 1:k) {
    if (i == k) {
      out <- paste0("fourd_",i,".nii.gz")
      fslmerge(subjList[(1 + (i-1)*break.subj):length.subj], direction="t", outfile=out, drop_dim=F)
      subMergeNames <- c(subMergeNames, out)
    } else {
      out <- paste0("fourd_",i,".nii.gz")
      fslmerge(subjList[(1 + (i-1)*break.subj):((i)*break.subj)], direction="t", outfile=out, drop_dim=F)
      subMergeNames <- c(subMergeNames, out)
    }
  }
  
  subMergeNames <- subMergeNames[-1]
  fslmerge(subMergeNames, direction="t", outfile="fourd.nii.gz")
  
  
  system('rm -f fourd_*.nii.gz')
  
  if (smooth > 0) {
    fslsmooth("fourd.nii.gz", sigma = smooth, outfile="fourd.nii.gz")
  } else {
    print("No smoothing done")
  }
  
  
} else {
  print("Skipping fourd image creation; Script will looking for file names fourd.nii.gz under first level directory")
}


system(paste0("scp ", maskName," ",OutDir, "/mask.nii.gz"), wait=T)
print("mask succesfully copied")


##############################################################################
################            Output summary files               ###############
##############################################################################

write.csv(subjData[, namePaths], paste0(namePaths,".csv"), row.names = F)
write.csv(subjData[, subjID], paste0(subjID,".csv"), row.names = F)

print("Succesfully wrote paths and id files")

##############################################################################
################        Make Output Directory                  ###############
##############################################################################


print("Creating output directory")
outName <- gsub("~", "", covsFormula)
outName <- gsub(" ", "", outName)
outName <- gsub("\\+","_",outName)
outName <- gsub("\\(","",outName)
outName <- gsub("\\)","",outName)
outName <- gsub(",","",outName)
outName <- gsub("\\.","",outName)
outName <- gsub("=","",outName)
outName <- gsub("\\*","and",outName)
outName <- gsub(":","and",outName)

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


mask<-readNIfTI(maskName)
imageIn<-readNIfTI(imageName)

###Time Series to Matrix Using oro.nifti
ts2matrix <- function(image, mask) {
  
  label <- sort(as.numeric(unique(matrix(mask@.Data))))
  
  if (length(label) == 2 && label[1] == 0 && label[2] == 1) {
    if (length(dim(image@.Data)) == 3 | dim(image@.Data)[4] == 1) {
      vector <- image@.Data[mask@.Data == 1]
      gc()
      return(vector)
      
    } else {
      temp <- matrix(image@.Data)[mask@.Data == 1]
      dim(temp) <- c(sum(mask@.Data), dim(image@.Data)[length(dim(image@.Data))])
      temp <- t(temp)
      
      temp <- as.data.frame(temp)
      names <- base::lapply(1:dim(temp)[2], function(x) { return(paste0("voxel",x))})
      names(temp) <- names
      gc()
      return(temp)
    }
    
  } else {
    gc()
    stop("Mask Image is not Binary")
  }
}



imageMat<-ts2matrix(imageIn,mask)

print("Fourd image and mask has been loaded")

##############################################################################
################           Preallocate output                  ###############
##############################################################################
pOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)
pAdjustedOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)
tOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)
zOut<-matrix(NA,nrow=dim(imageMat)[2],ncol=1)

print("Preallocate output done")

##############################################################################
################        Run model using parallel               ###############
##############################################################################
timeOn<-proc.time()

k <- 0
rm(k)
length.voxel <- ceiling(dim(imageMat)[2] / splits)

# We create a list of formulas for each voxel in our data. 
# Each element in the list will have formula with a different voxel as the dependent variable
print("Running Test Model")


m <- mclapply(1:10, function(x) {as.formula(paste(paste0("imageMat[,",x,"]"), covsFormula, sep=""))}, mc.cores = ncores)
test <- gam(formula = m[[1]], data=subjData, method="REML")


model <- mclapply(m, function(x) {
  foo <- summary(gam(formula = x, data=subjData, method="REML"))
  return(rbind(foo$p.table,foo$s.table))
}, mc.cores = ncores)


print("Test Models Done; Running Parallel Models")
for (k in 1:(splits)) {
  if (k == splits) {
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
      zOut[i,1]<- qnorm((model[[i]][which(rownames(model[[i]]) == variable),4] / 2), lower.tail=F)
    }
    
    pOutImage<-mask
    zOutImage<-mask
    
    pOutImage@.Data[mask==1@.Data]<-pOut
    zOutImage@.Data[mask==1@.Data]<-zOut
    
    pAdjustedOutImage<-mask
    pAdjustedOut <- stats::p.adjust(pOut, method=pAdjustMethod)
    pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
    
    var <- gsub("\\(", "", variable)
    var <- gsub("\\)", "", var)
    var <- gsub(",", "", var)
    var <- gsub("=", "", var)
    
    
    writeNIfTI(pOutImage,paste0("gamP_",var))
    writeNIfTI(zOutImage,paste0("gamZ_",var))
    if (pAdjustMethod != "none") {
      writeNIfTI(pAdjustedOutImage,paste0("gamPadjusted_",pAdjustMethod, "_",var)) 
    }
  }
  else {
    for(i in 1:length(model)){
      pOut[i,1]<-model[[i]][which(rownames(model[[i]]) == variable),4]
      zOut[i,1]<-sign(model[[i]][which(rownames(model[[i]]) == variable),3])*qnorm((model[[i]][which(rownames(model[[i]]) == variable),4] / 2), lower.tail=F)
      tOut[i,1]<-model[[i]][which(rownames(model[[i]]) == variable),3]
    }
    
    pOutImage<-mask
    zOutImage<-mask
    tOutImage<-mask
    
    pOutImage@.Data[mask==1@.Data]<-pOut
    zOutImage@.Data[mask==1@.Data]<-zOut
    tOutImage@.Data[mask==1@.Data]<-tOut
    
    pAdjustedOutImage<-mask
    pAdjustedOut <- stats::p.adjust(pOut, method=pAdjustMethod)
    pAdjustedOutImage@.Data[mask==1@.Data]<-pAdjustedOut
    
    var <- gsub("\\(", "", variable)
    var <- gsub("\\)", "", var)
    var <- gsub("\\*","and",var)
    var <- gsub(":","and",var)
    
    writeNIfTI(pOutImage,paste0("gamP_",var))
    writeNIfTI(zOutImage,paste0("gamZ_",var))
    writeNIfTI(tOutImage,paste0("gamT_",var))
    if (pAdjustMethod != "none") {
      writeNIfTI(pAdjustedOutImage,paste0("gamPadjusted_",pAdjustMethod, "_",var)) 
    }
    
  }
}

print("Write t-maps and p-maps is done")
