# Loaded data have gene expressions identified by EntrezID. As these are numeric, 
# some R functions have problem with this. To avoid this problem and also help 
# us with further analysis, such as conversion to gene set scores, we convert
# EntrezIDs to Gene Symbol.
# This function could also be used for conversion between different formats, if  
# provided data is not using Entrez ID

# BiocManager::install("org.Hs.eg.db")

# Load package for converting formats

library('org.Hs.eg.db') # Library for different gene namings conversion https://bioconductor.org/packages/release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf


#' Convert column names from EntrezID to Gene Symbols
#' @param inputData data with samples in rows and genes in columns for conversion
#' @return data with renamed columns
convertEntrezToGeneSymbol <- function(inputData){
  # Delete leading "X" present in some colnames 
  colnames(inputData) <- sub("X","", colnames(inputData), ignore.case = TRUE)
  
  # Set up map for conversion of entrezID to Gene Symbol
  geneSymbolsMap <- select(org.Hs.eg.db, keys = colnames(inputData), 
                           column = c("SYMBOL","ENTREZID"), keytype = "ENTREZID")
  # Convert Column names
  for (i in 1:dim(inputData)[2]){
    colnames(inputData)[i] <- 
      geneSymbolsMap$SYMBOL[colnames(inputData)[i]==geneSymbolsMap$ENTREZID]
  }
  return(inputData)
}
