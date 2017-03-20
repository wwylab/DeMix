#DeMix.Filter
############################################################

DeMix.Filter <- function(newt, groupid, zero_filter=TRUE, conc=0.8, fc=1.2)
{

  
  ## Filtering step 1.
  ## Input data are previously normalized at the data processing stage
  ## Don't use genes with count 0.
  genes.withoutzero=(apply(newt, 1, min)>0)
  
  if(zero_filter==TRUE){
    RNA=newt[genes.withoutzero==TRUE,]}
  
  
 # else {
 #   RNA <- newt}
  
  ## Filtering step 2
  ## Identify genes that satisfy the linearity Assumption with more than conc% probability for both up-regulated or down-regulated genes in tumor.
  
  Nsam=RNA[, (groupid==0)]
  Tsam=RNA[, (groupid==1)]
  
  Nmean <- apply(Nsam, 1, mean) ##mean gene expression for normals
  
  
  Cons <- rep(0, length(Nmean)) #logical vector for genes that are highly expressed in tumor
  Nega <- rep(0, length(Nmean)) #logical vector for genes that are highly expressed in normal
  
  for(i in 1:length(Nmean))
  {
    Cons[i] <- sum(Tsam[i,]>Nmean[i])
    Nega[i] <- sum(Tsam[i,]<Nmean[i])
  }
  cutoff_con <- round(sum(groupid==1)*conc)
  
  
  ## Filtering step 3
  ## Default fold-change is set at 1.2 and 1/1.2
  ## fc values need to be adjusted to make final dataset have around 2000~3000.
  
  efc=apply(Tsam, 1, mean)/(apply(Nsam, 1, mean)) ##ratio of the mean expression tumor vs. normal
  
  secfil= (apply(RNA, 1, max)<100000 & (efc> fc|efc< 1/fc) & ((Nega>cutoff_con) | (Cons>cutoff_con)) )
  
  inputdata <- RNA[secfil,]
  
  inputdata}