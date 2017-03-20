#DeMix.Normalization 
############################################################
DeMix.Normalization <- function(input,
                                NormalizationMethod="quantile",
                                groupid, ...) {
  row.id <- row.names(input)
  if(length(row.id) == 0) row.id = 1:nrow(input)
  newt <- as.matrix(input)
 
  ##source("https://bioconductor.org/biocLite.R")
  ##biocLite("DSS")
  require("DSS")
  
  colnames(newt) <- NULL
  rownames(newt) <- NULL
    
  # 0 denotes normal cell / 1 denotes tumor cell
  seqData <- newSeqCountSet(newt, groupid)
  
  ## Scale Normalization using 75th quantile
  ## Alternatively users can choose median or total for the scale normalization. Use "total" or "median" for the input of estNormFactors().
   seqData <- estNormFactors(seqData, NormalizationMethod) 
   k3 <- seqData@normalizationFactor
   mk3 <- median(k3)
   k3 <- k3/mk3

  # Calculating normalization factor
   for(i in 1:ncol(newt) ){
     newt[,i] = newt[,i]/k3[i]
    }
  
   row.names(newt) <- row.id
   return(newt)
}