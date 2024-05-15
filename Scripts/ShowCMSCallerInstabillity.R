# This script is used to generate a visualization of the instability of CMSCaller
# by running the classifier once to get classes, then discarding a selected class
# and rerunning the classifier on this reduced subset. The discarded class emerges 
# again.

# Load required libraries
library(CMScaller)
library(ggplot2)

# Load data
formatted_crc_data <- read.delim("formatted_crc_data.txt")
colnames(formatted_crc_data) <- sub("X","", colnames(formatted_crc_data), ignore.case = TRUE)

# Classify samples
res <- CMScaller(t(formatted_crc_data), doPlot=TRUE)
# CMS1 CMS2 CMS3 CMS4 <NA> 
# 432  859  508  987  446

# Classify again, this time ommit samples classified as CMS4 
dataWithoutCMS4 <- formatted_crc_data[is.na(res$prediction) | res$prediction != 'CMS4',]
res2 <- CMScaller(t(dataWithoutCMS4), doPlot=TRUE)
# CMS1 CMS2 CMS3 CMS4 <NA> 
# 312  555  348  677  353

# Create a data frame with classification count for plotting
predictionTable1 <- table(res$prediction)
predictionTable2 <- table(res2$prediction)

dataForPlot <- data.frame(
  Category = c(1,1,1,1,2,2,2,2), # Run of classifier
  CMS = c('CMS1','CMS2','CMS3','CMS4','CMS1','CMS2','CMS3','CMS4'),
  Values =  c(predictionTable1[[1]],predictionTable1[[2]],predictionTable1[[3]],predictionTable1[[4]],
              predictionTable2[[1]],predictionTable2[[2]],predictionTable2[[3]],predictionTable2[[4]])
)

customColors <- c("CMS1" = "orange", "CMS2" = "lightblue", "CMS3" = "purple", "CMS4" = "darkgreen")

ggplot(dataForPlot, aes(x = factor(Category), y = Values, fill = CMS)) +
  geom_bar(stat = "identity", width = 0.4) +
  scale_fill_manual(values = customColors) + 
  labs(title = "Two classifications by CMSCaller", x = "Classification", y = "Number of samples") +
  theme(
    text = element_text(size = 15),  # Adjust font size for all text elements
    axis.title = element_text(size = 15),  # Adjust font size for axis titles
    axis.text = element_text(size = 13),  # Adjust font size for axis labels
    legend.title = element_text(size = 15),  # Adjust font size for legend title
    legend.text = element_text(size = 13)  # Adjust font size for legend labels
  )
