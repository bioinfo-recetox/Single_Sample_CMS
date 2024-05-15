# This scripts contains the functions necessary to go from data and classification
# to a template usable with NTP. The steps used: 
# Match classifications to samples, as they are not in the same order and sometimes have a slightly 
# altered name (using "_" instead of " " and so on)

# Create a template based on RFE: 
# https://towardsdatascience.com/effective-feature-selection-recursive-feature-elimination-using-r-148ff998e4f7

#install.packages("caret")
#setwd("D:/Dokumenty/MUNI/BcPrace/Data")

library(caret) # for rfe
library(CMScaller) # for ntp
#library(CMSclassifier)

#head(emat) 
#head(res)



#' Order classifications the same order as samples by matching their sample names
#' @param inputData Data with samples in columns
#' @param inputClassifications Ground truth classifications containing (but not limited) to samples in inputData
matchClassificationsToSamples <- function(inputData, inputClassifications){
  
  # Samples in classifications have . replaced with -
  inputClassifications$sample <- gsub( "-","\\.", inputClassifications$sample)

  orderedClassification <- inputClassifications$CMS_final_network_plus_RFclassifier_in_nonconsensus_samples[match(colnames(inputData), inputClassifications$sample)]
  
  return(orderedClassification)
}

#' Create template from selected inputData and classifications using rfe.
#' @param inputData Data for template creation. Columns contain samples, rows contain genesets or gene expression scores.
#' @param orderedClassification Ground truth classification ordered the same way as inputData
#' @param rfeNumIter Number of iterations done in each step of RFE (as cross-validation)
#' @param rfeSizes Vector or number of features to try to reduce the number of features to.
#' @param printRFEResults Bool, print intermediary results
#' @param selectionThreshold Number>0, threshold importance to select a given feature. 
#' @return Template in the format of row with two columns - variable name and the class its important for
createTemplate <- function(inputData, orderedClassification, rfeNumIter=20, rfeSizes=c(5,15,25), 
                           printRFEResults=FALSE, selectionThreshold=0.1){
  
  ctrl <- rfeControl(functions = lmFuncs, allowParallel = T,
                     verbose=T, number = rfeNumIter)
  # ctrl <- rfeControl(functions = treebagFuncs, allowParallel = FALSE)
  
  print("Running RFE for feature selection (0 %)")
  CMS1Mask <- as.numeric(t(orderedClassification=="CMS1"))
  rfeResult1<-rfe(x = t(inputData), 
                 y = CMS1Mask, 
                 sizes = rfeSizes,
                 rfeControl = ctrl)
  
  print("RFE for CMS1 complete (25 %)")

  CMS2Mask <- as.numeric(t(orderedClassification=="CMS2"))
  rfeResult2<-rfe(x = t(inputData), 
                  y = CMS2Mask, 
                  sizes = rfeSizes,
                  rfeControl = ctrl)
  
  print("RFE for CMS2 complete (50 %)")
  
  CMS3Mask <- as.numeric(t(orderedClassification=="CMS3"))
  rfeResult3<-rfe(x = t(inputData), 
                  y = CMS3Mask, 
                  sizes = rfeSizes,
                  rfeControl = ctrl)
  
  print("RFE for CMS3 complete (75 %)")
  
  CMS4Mask <- as.numeric(t(orderedClassification=="CMS4"))
  rfeResult4<-rfe(x = t(inputData), 
                  y = CMS4Mask, 
                  sizes = rfeSizes,
                  rfeControl = ctrl)
  
  print("RFE for CMS4 complete (100 %)")
  
  print(paste("Number of selected predictors: CMS1:" ,length(predictors(rfeResult1)),
        "CMS2:", length(predictors(rfeResult2)),
        "CMS3:", length(predictors(rfeResult3)),
        "CMS4:", length(predictors(rfeResult4))))
  
  if(printRFEResults){
    print(rfeResult1)
    print(rfeResult2)
    print(rfeResult3)
    print(rfeResult4)
    # predictors(rfeResult1) # Shows RFE chosen predictors
    # ggplot(data = rfeResult1, metric = "RMSE") + theme_bw() # plot RMSE based on variables
  }
 
  print("Selecting features:")
  # For each row and each class compute means 
  CMS1Means <- rowMeans(inputData[,CMS1Mask==1]) 
  CMS2Means <- rowMeans(inputData[,CMS2Mask==1]) 
  CMS3Means <- rowMeans(inputData[,CMS3Mask==1]) 
  CMS4Means <- rowMeans(inputData[,CMS4Mask==1]) 
  CMSMeans <- cbind(CMS1Means,CMS2Means,CMS3Means,CMS4Means)
  
  selectedFeatures <- selectFeatures(rfeResult1,rfeResult2,rfeResult3,rfeResult4, importanceThreshold = selectionThreshold, CMSMeans)
  
  # Now finally create template, a df with 2 columns: probe and class
  newTemplate=data.frame(
    probe=selectedFeatures[1,],
    class=paste0("CMS",selectedFeatures[2,])
    )
  print("Template Created")
  return(newTemplate)
  
}

#' Algorithm to select features from RFE for the creation of a template
#' @param rfeResult1,rfeResult2,rfeResult3,rfeResult4 Results from the caret rfe function
#' @param importanceThreshold Threshold for comparison of feature importance
#' @param CMSMeans Table with rows containing variable names and cols containing means across CMS groups
#' @return A table with selected features and respective class
selectFeatures <- function(rfeResult1, rfeResult2, rfeResult3, rfeResult4, importanceThreshold=0.1, CMSMeans){
  
  # TODO: Do we always have all the features, or just the selected ones? In that case, a union would be needed
  #       and not present expressions would need to get 0 assigned.
  
  # Take into account feature importance: https://topepo.github.io/caret/recursive-feature-elimination.html  20.5.4
  
  # Take just variable + coeficient
  coefs1 <- rfeResult1$fit$coefficients
  coefs2 <- rfeResult2$fit$coefficients
  coefs3 <- rfeResult3$fit$coefficients
  coefs4 <- rfeResult4$fit$coefficients

  #effects1 <- rfeResult1$fit$residuals
  #effects2 <- rfeResult2$fit$residuals
  #effects3 <- rfeResult3$fit$residuals
  #effects4 <- rfeResult4$fit$residuals
  
  
  allFeatures <- unique(c(names(coefs1), names(coefs2), names(coefs3), names(coefs4)))
  
  templateSets<-c()
  templateClassifications<-c()

  # Loop over every gene set name and select to which group it belongs
  # based on maximum coefficient 
  
  for(i in 2:length(allFeatures)){   # i=1 is intercept
    setName <- allFeatures[i]
    val1 <- 0
    val2 <- 0
    val3 <- 0
    val4 <- 0
    # Feature doesn't have to be in a set, prevent subscript out of errors 
    try({val1 <- coefs1[[setName]]},silent=T)
    try({val2 <- coefs2[[setName]]},silent=T)
    try({val3 <- coefs3[[setName]]},silent=T)
    try({val4 <- coefs4[[setName]]},silent=T)
    # Features that are negatively expressed are multiplied by 0
    val1 <- val1 *( CMSMeans[setName,"CMS1Means"] > 0)
    val2 <- val2 *( CMSMeans[setName,"CMS2Means"] > 0)
    val3 <- val3 *( CMSMeans[setName,"CMS3Means"] > 0)
    val4 <- val4 *( CMSMeans[setName,"CMS4Means"] > 0)
    
    # If the feature is not important in any set, skip it
    if(max(c(val1,val2,val3,val4)) < importanceThreshold){
      next  
    }
    # select in which set this feature is the most important
    highestSet <- which.max(c(val1,val2,val3,val4))
    templateSets<-c(templateSets, setName)
    templateClassifications<-c(templateClassifications, highestSet)
  }
  
  return(rbind(templateSets,templateClassifications))
} 
