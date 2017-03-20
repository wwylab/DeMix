\name{DeMix.Normalization}
\alias{DeMix.Normalization}
\title{
   Normalization for DeMix
}
\usage{
DeMix.Normalization(input,NormalizationMethod, groupid)
}
\description{
A normalization function for quantile/total/median normalization.
The preparation work -- quantile normalization before DeMix (Bayesian deconvolution models for both RNA-seq read counts and microarray expressions).} 

\arguments{
\item{input}{
yourDataInput: A matrix of gene expressions before normalization and selection. The GxS matrix where G is the
number of genes and S is the number of samples.}

\item{NormalizationMethod}{Method for normalization. Here we use scale normalization. The value used to rescale can be the 75th quantile, i.e. "quantile" which is the default. The median, i.e. "median", or the column sum, i.e. "total" can also be used as alternative methods for the scale normalization. }


\item{groupid}{A vector of indicators of historical normal samples and mixed tumor samples. We use 1 to denote the mixed tumor samples and 0 to denote the normal samples.}

}

\value{
 \item{newtfile}{A matrix of gene expression levels. The GxS matrix where G is the number of genes and S is the number of samples.}}

\author{
Jaeil Ahn
}
\references{
Ahn J, Yuan Y, Parmigiani G, Suraokar MB, Diao L, Wistuba II, and Wang W. Bioinformatics 2013 doi: 10.1093/bioinformatics/btt301.

Wu H, Wang C and Wu Z (2013). “A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data.” Biostatistics. doi: 10.1093/biostatistics/kxs033.
}
\keyword{Normalization}
\seealso{
 http://bioinformatics.mdanderson.org/main/DeMix:Overview
 }
\examples{
 ##perform normalization on a simulated dataset
 ##there are 10 normal samples and 10 tumor samples 	
 data(simulS1)
 input <- simulS1
 newtfile <- DeMix.Normalization(input, "quantile",c(rep(0, 10), rep(1, 10)))
}
