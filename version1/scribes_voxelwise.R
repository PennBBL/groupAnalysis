##############################################################################
################                 Option List                   ###############
##############################################################################

#library Oronifti
library(oro.nifti)

maskName <- paste0(OutDir,"/mask.nii.gz")
imageName <- paste0(OutDir,"/fourd.nii.gz")

mask<-readNIfTI(maskName)
imageIn<-readNIfTI(imageName)

tsExample <- ts2matrix(imageIn, mask)
test <- imageIn@.Data[,,,1]
which(test == -7.255098338126042101237)

mask2<-antsImageRead(maskName,3)
imageIn2<-antsImageRead(imageName,4)
imageMat<-timeseries2matrix(imageIn,mask)

test2 <- imageIn2[1:91,1:109,1:91,1]
which(test2 == -7.2550983428955078125)

image1 <- "/import/monstrum/eons_xnat/subjects/80010_2894/7_bbl1_frac2back1_231/stats/002894_bbl1_frac2back1_231_SEQ07_eons_frac_nobehave_stats_RUN01.feat/reg_std_ants/cope4_2mm.nii.gz"

image <- readNIfTI2(image1)
image2 <- antsImageRead(image1, 4)
which(image@.Data == -7.255098338126042101237)
which(image2 == -7.2550983428955078125)

imageMat <- tsmeanCluster(imageIn, mask)



##############################################################################
################                 Option List                   ###############
##############################################################################


gamPtable <- function(imageMatrix, covsFormula, subjData) {
  
  foo <- summary(gam(formula = covsFormula, data=subjData, method="REML"))
  return(rbind(foo$p.table,foo$s.table))
  
}


gammPtable <- function(imageMatrix, covsFormula, subjData, randomFormula) {
  
  foo <- summary(gamm4(formula = covsFormula, data=subjData, REML=T, random = as.formula(randomFormula))$gam)
  return(rbind(foo$p.table,foo$s.table))
  
}

outPutMaps <- function(model) {
  for (j in 1:dim(model[[1]])[1]) {
    
    pOut<-matrix(NA,nrow=length(model),ncol=1)
    tOut<-matrix(NA,nrow=length(model),ncol=1)
    zOut<-matrix(NA,nrow=length(model),ncol=1)
    
    variable <- rownames(model[[1]])[j]
    
    if (grepl("s(", variable, fixed=T)) {
      
      for(i in 1:length(model)){
        pOut[i,1]<- model[[i]][which(rownames(model[[i]]) == variable),4]
        zOut[i,1]<- qnorm(model[[i]][which(rownames(model[[i]]) == variable),4], lower.tail=F)
      }
      
      pOutImage<-mask
      pOutImage@.Data[mask@.Data==1]<-pOut
      
      zOutImage<-mask
      zOutImage@.Data[mask@.Data==1]<-zOut
      
      var <- gsub("\\(", "", variable)
      var <- gsub("\\)", "", var)
      var <- gsub(",", "", var)
      var <- gsub("=", "", var)
      
      
      
      writenii(pOutImage,paste0("gamm4P_",var,".nii.gz"))
      writenii(zOutImage,paste0("gamm4Z_",var,".nii.gz"))
    }
    else {
      for(i in 1:length(model)){
        pOut[i,1]<-model[[i]][which(rownames(model[[i]]) == variable),4]
        zOut[i,1]<-qnorm(model[[i]][which(rownames(model[[i]]) == variable),4], lower.tail=F)
        tOut[i,1]<-model[[i]][which(rownames(model[[i]]) == variable),3]
      }
      
      pOutImage<-mask
      zOutImage<-mask
      tOutImage<-mask
      pOutImage@.Data[mask==1@.Data]<-pOut
      zOutImage@.Data[mask==1@.Data]<-zOut
      tOutImage@.Data[mask==1@.Data]<-tOut
      
      var <- gsub("\\(", "", variable)
      var <- gsub("\\)", "", var)
      var <- gsub("\\*","and",var)
      var <- gsub(":","and",var)
      
      writenii(pOutImage,paste0("gamm4P_",var,".nii.gz"))
      writenii(zOutImage,paste0("gamm4Z_",var,".nii.gz"))
      writenii(tOutImage,paste0("gamm4T_",var,".nii.gz"))
      
    }
  }
}
