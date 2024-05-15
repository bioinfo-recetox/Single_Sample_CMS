# Convert gene expression to pathway activation scores. 
# Uses GSEA MSigDB -- Hallmark gene sets. This set contains 50 well
# described biological states or processes.
# Source: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp#H
# To these are added additional gene sets from literature, as described in the thesis.

# Uses the Bioconductor/GSVA package, specifically the single sample ss-GSVA 

# install GSVA package https://bioconductor.org/packages/release/bioc/manuals/GSVA/man/GSVA.pdf 
# BiocManager::install("GSVA")
# install.packages("msigdbr")

library(GSVA)
library(msigdbr)

# Link to set up: https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_03_gsva.html#4_Gene_set_variation_analysis_-_RNA-Seq

#' Convert matrix with gene expression data to a matrix with pathway activation scores 
#' 
#' @param inputData Matrix, where rows contain Gene symbols and cols sample IDs
#' @returns A converted matrix with rows reduced to hallmark pathways.
ConvertToPathways <- function(inputData){
  # Import hallmark dataset
  hallmarkGeneSets <- msigdbr::msigdbr(
    species = "Homo sapiens", # Can change this to what species you need
    category = "H" # Only hallmark gene sets
  )

  # length(unique(hallmarkGeneSets$gs_name)) # Hallmark DB should have 50 pathways.
  
  # Reformat hallmarkGeneSets into lists, using gene symbols.
  listHallmarkGeneSets <- split(hallmarkGeneSets$human_gene_symbol, hallmarkGeneSets$gs_name)
  
  # Load additional pathways from literature
  load("./Additional_gene_sets/pelka_genesets.rdata")
  
  # Join lists: 
  listGeneSets <- c(listHallmarkGeneSets, epithelial.cells.programs, immune.cells.programs, stromal.cells.programs)
  
  # Use single-sample function gsva::gsva to convert expression to gene set scores
  # using supplied list of gene sets
  parametersForGSVA <- gsvaParam(exprData = inputData, geneSets = listGeneSets)
  pathwayScores <- gsva(parametersForGSVA)
  
  
  # Remove special characters from rownames, which are not compatible with RFE
  rownames(pathwayScores) <- gsub( "\\-","_", rownames(pathwayScores))
  rownames(pathwayScores) <- gsub( "\\.","_", rownames(pathwayScores))
  rownames(pathwayScores) <- gsub( "\\,","_", rownames(pathwayScores))
  rownames(pathwayScores) <- gsub( "\\(","", rownames(pathwayScores))
  rownames(pathwayScores) <- gsub( "\\)","", rownames(pathwayScores))
  rownames(pathwayScores) <- gsub( "\\/","", rownames(pathwayScores))
  
  # Convert pathway scores table to numeric
  originalRowNames <- row.names(pathwayScores)
  pathwayScores <- apply(pathwayScores,2, as.numeric)
  rownames(pathwayScores) <- c(originalRowNames)
  
  return(pathwayScores);
}

