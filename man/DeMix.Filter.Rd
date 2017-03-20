\name{DeMix.Filter}
\alias{DeMix.Filter}
\docType{package}
\title{
    Pre-selection of a subset of genes for the estimation of mixing proportions
}

\description{
To avoid unnecessary computational burden from analyzing all 20,000 genes as well as added noise from genes with similar expression levels in the N and T components, gene pre-selection is a critical step for the accurate estimation of mixing proportions.
}

\usage{
DeMix.Filter(data=newt, groupid, zero_filter=TRUE, conc=0.8, fc=1.2)
}


\arguments{
\item{data}{
data: A matrix of normalized gene expressions. The GxS matrix where G is the
number of genes and S is the number of samples.}

\item{groupid}{A vector of tissue component, 0 for normal and 1 for tumor samples}

\item{zero_filter}{An indicator of whether the zero counts will be removed. The default is True, which indicates remove.}


\item{conc}{A cut-off value used for checking the linearity assumption. More than 0.8 is recommended.
 }

\item{fc}{A scalar value to provide a fold-change cut-off, in order to identify differentially expressed genes to be used in proportion estimation. The default is 1.2.}



}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{newtfile}{A new matrix of gene expressions after the gene pre-selection. The G'xS matrix where G' is the number of remained genes and S is the number of samples.}

}

\author{
Jaeil Ahn
}
\references{
Ahn J, Yuan Y, Parmigiani G, Suraokar MB, Diao L, Wistuba II, and Wang W. Bioinformatics 2013 doi: 10.1093/bioinformatics/btt301.

Wu H, Wang C and Wu Z (2013). “A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data.” Biostatistics. doi: 10.1093/biostatistics/kxs033.
}
\keyword{Filter}
\seealso{
 http://bioinformatics.mdanderson.org/main/DeMix:Overview
 }


\examples{
data(simulS1)
input <- as.matrix(simulS1)
newtfile <- DeMix::DeMix.Filter(input,c(rep(0, 10), rep(1, 10)), zero_filter=TRUE, conc=0.8, fc=1.2)
}




