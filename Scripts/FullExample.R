# Include needed packages
library(ggplot2)
library(CMScaller)
library(org.Hs.eg.db)

# Source needed functions from external files 
source("CreateTemplate.r")
source("GenesToPathways.r")
source("ConvertFormats.r")
source("FilterExpressionData.R")

# Download from synapse requires a login
# synapseID <- "syn4983432"
# synLogin(email= ... , authToken= ...)
# formattedCrcData <- downloadSynapseData(synapseID ,saveLocally = F)

# (15. 05. 2024) This throws error: File syn4983432 not found in Synapse
# The file can be downloaded manually from https://www.synapse.org/#!Synapse:syn4983432 

# Import pre-saved data, as download takes long
getwd()
formattedCrcData <- read.delim("formatted_crc_data.txt")

# Convert EntrezID to GeneSymbol
formattedCrcData <- convertEntrezToGeneSymbol(formattedCrcData)

# Convert expression data to pathway activation scores (takes about 5 minuts)
pathwaysFromFormatedData <- ConvertToPathways(t(formattedCrcData))

# save pathway scores
write.table(pathwaysFromFormatedData,file="pathwayScoresAll.table")

# load pathway scores
# Scores with just Hallmark
# pathwaysFromFormatedData2 <- read.table("pathwaysFromFormatedData.table")
# Added Epithelilal gene sets
# pathwaysFromFormatedData2 <- read.table("pathwayScoresWithAddedEpithelial.table")
# pathwaysFromFormatedData2 <- data.matrix(pathwaysFromFormatedData2)
# Added all gene sets
pathwaysFromFormatedData2 <- read.table("pathwayScoresAll.table")
pathwaysFromFormatedData2 <- data.matrix(pathwaysFromFormatedData2)


# Download / Load original CRCSC classifications 
# synapseID <- "syn4978511"
# originalClassifications <- downloadSynapseData(synapseID ,saveLocally = F)
originalClassifications <- read.delim("cms_labels_public_all.txt")

# Match original classifications to samples 
matchedClassifications <- matchClassificationsToSamples(pathwaysFromFormatedData2, originalClassifications)
matchedClassificationsForExpression <- matchClassificationsToSamples(t(formattedCrcData), originalClassifications)


# Separate training and test data.
set.seed(201)
inTrain <- createDataPartition(matchedClassifications, p = .80, list = FALSE)[,1]
inTrainForExpression <- createDataPartition(matchedClassificationsForExpression, p = .50, list = FALSE)[,1]

# Now create template for Gene Sets (time-intensive, especially with high iteration counts. Cached results can be loaded below)
geneSetsTemplate <- createTemplate(pathwaysFromFormatedData2[,inTrain], matchedClassifications[inTrain],
                               rfeNumIter = 50, rfeSizes = c(15,35,45,55), printRFEResults = F, selectionThreshold=0.05)

# We can view the template and number of features for each subclass
geneSetsTemplate
table(geneSetsTemplate$class)

save(geneSetsTemplate, file = "./results/allGeneSetsTemplate.rdata")
# load("./results/allGeneSetsTemplate.rdata")

# And create template for Gene expressions
if(FALSE){ # Takes a long time to compute
exprTemplate <- createTemplate((formattedCrcData)[,inTrainForExpression], matchedClassificationsForExpression[inTrainForExpression],
                              rfeNumIter = 20, rfeSizes = c(50), printRFEResults = TRUE, maxParams = 50)

}

# Use created template for classification of test and all variables

newPredictions <- ntp(pathwaysFromFormatedData2[,-inTrain], geneSetsTemplate, doPlot=TRUE, nPerm=100, distance = "cosine")
newPredictionsAll <- ntp(pathwaysFromFormatedData2, geneSetsTemplate, doPlot=TRUE, nPerm=100, distance = "cosine")

#pValueThreshold <- 0.05
newPredictions$prediction <- as.character(newPredictions$prediction)
#newPredictions$prediction <- ifelse(newPredictions$p.value > pValueThreshold, "NOLBL", newPredictions$prediction)

# And compare with original predictions in train data 

totalAccuracy <- sum(newPredictions$prediction == matchedClassifications[-inTrain])/ length(matchedClassifications[-inTrain])
totalAccuracy 
# 0.7065217, 0.6599379,  0.6490683

# Create a confusion matrix, which also gives a confidence interval for the accuracy
confusionMatrix(as.factor(newPredictions$prediction), as.factor(matchedClassifications[-inTrain]))

# We can compare the accuracy to the accuracy of CMSCaller measured the same way (% of identity) 

CMSCallerRes <- CMScaller(t(formattedCrcData[-inTrain,]),rowNames="symbol")
#replace NA values
CMSCallerRes$prediction <- ifelse(is.na(CMSCallerRes$prediction), "NOLBL", CMSCallerRes$prediction)
CMSCallerRes$prediction <- paste0("CMS",CMSCallerRes$prediction)
CMSCallerRes$prediction <- ifelse((CMSCallerRes$prediction == "CMSNOLBL"), "NOLBL", CMSCallerRes$prediction)

CMSCallerAccuracy <- sum(CMSCallerRes$prediction == matchedClassifications[-inTrain])/ length(matchedClassifications[-inTrain])
CMSCallerAccuracy 
# 0.7236025

# Create a confusion matrix, which also gives a confidence interval for the accuracy
confusionMatrix(as.factor(CMSCallerRes$prediction), as.factor(matchedClassifications[-inTrain]))



# ----------- Gene expressions ----------------
# Data contains too many variables for RFE
# First, we need to filter them by SD: 
ratioToDiscard <- 0.7
filteredCrcData <- removeLowSDData(t(formattedCrcData), ratioToDiscard)

# Convert to Gene Symbol, as numeric rownames are a problem, remove '-' in row names
# filteredCrcData <- t(convertEntrezToGeneSymbol(t(filteredCrcData)))
rownames(filteredCrcData) <- gsub("-", "_", rownames(filteredCrcData))
# For this, we need to normalize the data
summary(t(filteredCrcData)[,1:10])
scaledData <- t(scale(t(filteredCrcData)))
summary(t(scaledData)[,1:10])

geneExpressionTemplate <- createTemplate(scaledData[,inTrain], matchedClassifications[inTrain],
                                   rfeNumIter = 10, rfeSizes = c(60,100), printRFEResults = F, selectionThreshold=0.015)
table(geneExpressionTemplate$class)
# save(geneExpressionTemplate, file = "./results/allGeneExpressionTemplate.rdata")
# load("./results/allGeneExpressionTemplate.rdata")

# Use created template for classification of test and all variables
newPredictions2 <- ntp(scaledData[,-inTrain], geneExpressionTemplate, doPlot=TRUE, nPerm=100, distance = "cosine")
newPredictionsAll <- ntp(scaledData, geneExpressionTemplate, doPlot=TRUE, nPerm=100, distance = "cosine")

totalAccuracy2 <- sum(newPredictions2$prediction == matchedClassifications[-inTrain])/ length(matchedClassifications[-inTrain])
totalAccuracy2 
#  0.7236025, 0.6956522

# And again, compute confusion matrix and CI
confusionMatrix(as.factor(newPredictions2$prediction), as.factor(matchedClassifications[-inTrain]))
