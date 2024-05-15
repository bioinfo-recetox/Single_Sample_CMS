# This file is used to reduce the number of genes 

# Not run with sourcing file
if(FALSE){
  
formatted_crc_data <- read.delim("formatted_crc_data.txt")

# medians for visualisation

medians <- apply(formatted_crc_data, 1, median)
hist(medians)

# Compute SD across genes
standardDeviations <-  apply(formatted_crc_data, 1, sd)

# Visualise SD
hist(standardDeviations, main="Standard deviations in expression data", 
     ylim=c(0,800), xlim=c(1.2,2),
     xlab="SD", ylab="Count")

# Compute the value at threshold quantile
thresholdQuantile <- 0.4 # 0.4 -> discard 40 % of data
thresholdSD <- quantile(standardDeviations, probs = thresholdQuantile)

# Take only where SD is over threshold
filteredCrcData <- formatted_crc_data[standardDeviations > thresholdSD,]
}


# Remove variables with lowest SDs from data. 
#
# @param data where rows contain variables to filter
# @param ratioToDiscard value between 0 and 1, what proportion of values gets discarded
removeLowSDData <- function(data, ratioToDiscard){
  # Compute SD across genes
  standardDeviations <-  apply(data, 1, sd)
  
  # Compute the value at threshold quantile
  thresholdQuantile <- ratioToDiscard 
  thresholdSD <- quantile(standardDeviations, probs = thresholdQuantile)
  
  # Take only where SD is over threshold SD
  filteredCrcData <- data[standardDeviations > thresholdSD,]
  
  return(filteredCrcData)
}
