\name{DeMix.Model}
\alias{DeMix.Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Decovolution for mixed cancer transcriptomes using raw measured data from RNAseq or microarray platforms
}
\description{
In this stage, the deconvolution is designed for a probeset. We focus on the estimation of \eqn{\pi} , the purity for the unknown tumor component. The following three steps comprise the selection procedure of this subset.
Step 1: The probes are selected with small standard deviations of the log2- transformed intensity for normal tissues from historical samples. Step 2: Those selected probes are sorted by their standard deviations of raw intensity and those
ones with the largest standard deviations are prepared for the next step. Step 3: The DeMix is run on all probes to estimate \eqn{\mu} and \eqn{\pi} for each probeset.
}
\usage{
DeMix.Model(inputdata, groupid, givenpi, NormalizationStep, nStage, niter)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputdata}{
A matrix of gene expressions from mixed tumor and normal samples. The GxS matrix where G is the
number of genes and S is the number of samples.}
  
  \item{groupid}{
A vector of indicators of historical normal and mixed tumor samples. We use 1 to denote mixed tumor samples and 0 to denote normal samples.}

  \item{givenpi}{An indicator or a vector. We use 1 to denote the \eqn{\pi} for each sample is unknown. If the proportion of tumor component is known, givenpi is the vector to denote the proportion.}
  
  \item{NormalizationStep}{An indicator. We use 1 to denote the need for normaliztion. We use 0 to denote no need for normalization. }
  
  \item{NormalizationMethod}{A charactor, if NormalizationStep==1, then the corresponding method used for scale normalization. The value used to rescale can be the 75th quantile, i.e. "quantile" which is the default. The median, i.e. "median", or the column sum, i.e. "total" can also be used as alternative methods for the scale normalization. }
 
  \item{nStage}{The number of nStage must be 1, 2 or 3. The default is 3. 
  	
  	If nStage is 1, DeMix will only perform the estimation of \eqn{\pi}, the proportion for the unknown component, i.e. tumor. DeMix first select a subset of genes/probes using \link{DeMix.filter}, in order to reduce the computational burden and increase accuracy in the proportion estimation. 
  	
  	The following three steps comprise the selection procedure of this subset.
Step 1: The probes are selected with small standard deviations of the log2- transformed intensity for normal tissues from historical samples. Step 2: Those selected probes are sorted by their standard deviations of raw intensity and those
ones with the largest standard deviations are prepared for the next step. Step 3: The DeMix is run in those pre-selected genes/probes to estimate \eqn{\mu} and \eqn{\pi}. 

  If nStage is 2, DeMix estimates individual expression profiles for both tumor and stromal tissues for a given gene subset (or whole gene set), using the known proportion of tumor component from givenpi input. 

	If nStage is 3, DeMix first runs stage 1 and then further estimates individual expression profiles for both tumor and stromal tissues for a given gene subset (or whole gene set), using the proportion estimates of tumor component from nStage 1. 
}
 
 \item{seeds}{Seeds set for the random number generation for MCMC sampler. The default for seeds is c(629555906, 921927245, 1265635378).
}

 \item{niter}{The iteration number of MCMC process. The default for niter is 1500.
}
 
  }

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{pi}{Estimated tumor proportion: 6 rows each sample: 3 rows from deconvolution by assuming Poisson distribution, 3 rows form deconvolution by assuming Negative Binomial distribution;
                 posterior means are given on the 1st and 4th row; 95\% Confidence Interval are given in the 2nd, 3rd, 5th, 6th row}
                 
  \item{decovMu}{Estimated \eqn{\mu} for Negative Binomial distribution,\eqn{\mu_{1}} is for normal samples and \eqn{\mu_{2}} is for tumor samples.}
  \item{decovExpr_normal}{A matrix with deconvolved expression profile of normal. For nStage=2, the matrix contains the full set of genes. For nStage=1, the matrix contains a pre-selected set of genes based on \link{DeMix.Filter}.}
   \item{decovExpr_tumor}{A matrix with deconvolved exprssion profile of tumor. For nStage=2,the matrix contains the full set of genes. For nStage=1, the matrix contains a pre-selected set of genes based on \link{DeMix.Filter}.}
  
   \item{post_pi}{A matrix with posterior samples for \eqn{\pi}. Only 10 designated posterior samples are stored, while the default is top 10 genes, from gene 1 to gene 10 by setting pgenes=c(1,2,3,4,5,6,7,8,9,10).}
   
   \item{post_decovMu}{A matrix with posterior samples for \eqn{\mu}. Only 10 designated posterior samples are stored, while the default is top 10 genes, from gene 1 to gene 10 by setting pgenes=c(1,2,3,4,5,6,7,8,9,10).}
   }
\note{
We propose a Bayesian deconvolution models (DeMix-Bayes) for both RNA-seq read counts and microarray expressions. Using an MCMC algorithm, DeMix-Bayes address two
challenges: 1) estimation of both tumor proportion and tumor-specific expression, when
neither is known a priori, 2) estimation of individualized expression profiles for both tumor
and stromal tissues. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\author{
Jaeil Ahn
}
\references{
Ahn J, Yuan Y, Parmigiani G, Suraokar MB, Diao L, Wistuba II, and Wang W. Bioinformatics 2013 doi: 10.1093/bioinformatics/btt301.

Wu H, Wang C and Wu Z (2013). “A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data.” Biostatistics. doi: 10.1093/biostatistics/kxs033.
}
\keyword{DeMix.Model}
\seealso{
 http://bioinformatics.mdanderson.org/main/DeMix:Overview
 }


\examples{
###################################################################################
## We focus on the estimation of the purity for the unknown tumor component,
## then further estimates individual expression profiles for both 
## tumor and stromal tissues for a given gene subset (or whole gene set), using the ## proportion estimates of tumor component.
data(simulS1)
input <- as.matrix(simulS1)
groupid <- c(rep(0, 10), rep(1, 10))
givenpi <- 1
NormalizationStep <- 1
NormalizationMethod <- "quantile"
nStage <- 3
res <- DeMix::DeMix.Model(input, groupid, givenpi, NormalizationStep, NormalizationMethod,nStage)
##print DeMix output and compare with the true proportions
##print(res$pi)
##print(res$post_pi)
##print(res$post_decovMu)
##print(Prop.sim)

###################################################################################
## The default iteration number of MCMC process is 1500. We also can set it as 1000 ## or other values to reduce the iteration time based on the example above.
##data(simulS1)
##input <- as.matrix(simulS1)
##groupid <- c(rep(0, 10), rep(1, 10))
##NormalizationStep <- 1
##NormalizationMethod <- "quantile"
##nStage <- 3
##niter <- 1000
##res <- DeMix::DeMix.Model(input, groupid, givenpi, NormalizationStep, Normalizatio##nMethod,nStage,niter)

###################################################################################
## If the uses already have the known proportion of tumor component from givenpi 
## input, we can set nStage as 2 in order to obtain DeMix estimates individual 
## expression profiles for both tumor and stromal tissues for a given gene subset. 
##data(simulS1)
##input <- as.matrix(simulS1)
##data(simu_givenpi)
##simu_givenpi <- as.matrix(simu_givenpi)
##groupid <- c(rep(0, 10), rep(1, 10))
##givenpi <- 2
##NormalizationStep <- 1
##NormalizationMethod <- "quantile"
##nStage <- 2
##niter <- 1000
##res <- DeMix::DeMix.Model(input, groupid, givenpi, NormalizationStep, Normalizatio##nMethod,nStage)
##print DeMix output and compare with the true proportions
##print(res$pi)
##print(res$post_pi)
##print(res$post_decovMu)
##print(Prop.sim)

###################################################################################
## When running nStage 2, if you have genes more than 1000, in order to save running ## time, we suggest to divide into small groups. As an example below, suppose we 
## have 3000 genes in input file, we can split the input firstly.
##input <- read.table("./yourInput.txt")
##input1 <- input[1:1000, ]
##input2 <- input[1001:2000, ]
##input3 <- input[2001:3000, ]
##data(simu_givenpi)
##simu_givenpi <- as.matrix(simu_givenpi)
##groupid <- c(rep(0, 10), rep(1, 10))
##givenpi <- 2
##NormalizationStep <- 1
##NormalizationMethod <- "quantile"
##nStage <- 2
##niter <- 1000
##res1 <- DeMix::DeMix.Model(input1, groupid, givenpi, NormalizationStep, Normalizat##ionMethod,nStage)
##res2 <- DeMix::DeMix.Model(input2, groupid, givenpi, NormalizationStep, Normalizat##ionMethod,nStage)
##res3 <- DeMix::DeMix.Model(input3, groupid, givenpi, NormalizationStep, Normalizat##ionMethod,nStage)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DeMix}% use one of  RShowDoc("KEYWORDS")
\keyword{DeMix}% __ONLY ONE__ keyword per line
