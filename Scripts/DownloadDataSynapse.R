# This file contains functions to download needed data from the Synapse website
# Data can be downloaded based on their synapse ID 


# Synapse project: https://www.synapse.org/#!Synapse:syn2623706/files/

library(synapser)

# synLogin requires email and authToken. Either provide them as rows of a textfile, add them straight into the synLogin
# or download data from the website and not through code
loginInfo <- readLines("Synapse_credentials.txt")
synLogin(email=loginInfo[1], authToken=loginInfo[2])
library("annotate")

# syn4978510 clinical_molecular_public_all.txt 
# syn4978511 cms_labels_public_all.txt
# syn2321865 curated genesets

# syn4983432 formatted_crc_data.txt


#' Function that returns downloaded data from Synapse based on synId. Optionally can save this data locally, 
#' so it does not need to be redownloaded again. 
#' @param synId Synapse ID of the file that should be downloaded
#' @param saveLocally Bool, whether the file should be saved into the working directory under the synId
#' @return data matrix
#' @note This function requires a Synapse login. This is done via the function synLogin(email, authToken)
downloadSynapseData <- function(synId, saveLocally=F){
  # Download
  tmp <- synGet(synId)
  t <- read.delim(tmp$path, as.is=T, header=F)
  rns <- t[-1, 1]
  cns <- t[1, -1]
  t <- t[-1, -1]
  rownames(t) <- rns
  colnames(t) <- cns

  if(saveLocally){
    write.table(t, file=synId)
  }
  return(t)
}
