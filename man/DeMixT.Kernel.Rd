\name{DeMixT.Kernel}
\alias{DeMixT.Kernel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Kernel function for optimizing parameters and hidden variables in DeMixT
}
\description{
This function is invoked by DeMixT.S1 and DeMixT.S2 to finish parameter estimation and expression deconvolution.
}
\usage{
DeMixT.Kernel(inputdata, groupid, nhavepi, givenpi, givenpiT, niter, ninteg)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inputdata}{
A matrix of expression data (e.g gene expressions) from reference (e.g. normal) and mixed samples (e.g. mixed tumor samples). It is a \eqn{GxS} matrix where \eqn{G} is the number of genes and \eqn{S} is the number of samples including reference and mixed samples. Samples with the same tissue type should be placed together in columns (e.g. cbind(normal samples, mixed tumor samples)}
  \item{groupid}{
A vector of indicators to denote if the corresponding samples are reference samples or mixed tumor samples. DeMixT is able to deconvolve mixed tumor samples with at most three components. We use 1 and 2 to denote the samples referencing the first and the second known component in mixed tumor samples. We use 3 to indicate mixed tumor samples prepared to be deconvolved. For example, in two-component deconvolution, we have c(1,1,...,3,3) and in three-component deconvolution, we have c(1,1,...,2,2,...,3,3).}
  \item{nhavepi}{If it is set to 0, then deconvolution is performed without any given proportions; if set to 1, deconvolution with given proportions for the first and the second known component is run; if set to 2, deconvolution is run with given tumor proportions. This option helps to do deconvolution in different settings. Because in estimation of component-specific proportions, we just use a subset of genes ; so when it is required to deconvolve another subset of genes, we just easily plug back our estimated proportions by setting this option to 1. In our two-step estimation strategy in a three-component setting, this option is set to 2 to implement the second step.}
 \item{givenpi}{\eqn{S_{T_N}}-Vector of proportions. Given the number of mixed tumor samples is \eqn{S_T(S_T<S)}, \eqn{S_{T_N}} is set to \eqn{2*S_T} in a three-component setting and \eqn{S_T} in a two-component setting. When nhavepi is 1, it is fixed with the given proportions for the first and the second known component of mixed tumor samples, or just for one known component when there is just one type of reference tissues. It has the form of Vector \eqn{\pi^1_{N_1}},\eqn{pi^2_{N_1}},...,\eqn{\pi^{S_T}_{N_1}},\eqn{\pi^1_{N_2}},\eqn{\pi^2_{N_2}},...,\eqn{\pi^{S_T}_{N_2}}.}
\item{givenpiT}{\eqn{S_T}-Vector of proportions. When nhavepi is set to 2, givenpiT is fixed with given proportions for unknown component of mixed tumor samples. This option is used when we adopt a two-step estimation strategy in deconvolution. It has the form of Vector \eqn{\pi^1_T}, \eqn{\pi^2_T}, \eqn{\cdots}, \eqn{\pi^{S_T}_T}. If option is not 2, this vector can be given with any element.}
\item{niter}{The number of iterations used in the algorithm of iterated conditional modes. A larger value can better guarantee the convergence in estimation but increase the computation time.}
\item{ninteg}{The number of bins used in numerical integration for computing complete likelihood. A larger value can increase accuracy in estimation but also increase the running time. Especially in three-component deconvolution, the increase of number of bins can greatly lengthen the running time.}

}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{pi}{Matrix of estimated proportions for each known component. The first row corresponds to the proportion estimate of each sample for the first known component (groupid = 1) and the second row corresponds to that for the second known component (groupid = 2)} 
 \item{decovExpr}{A matrix of deconvolved expression profiles corresponding to unknown (e.g tumor) component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample.}           
 \item{decovMu}{Estimated \eqn{\mu} of log2-normal distribution for tumor component.}
 \item{decovSigma}{Estimated \eqn{\sigma} of log2-normal distribution for tumor component}
 \item{pi1}{An \eqn{S_T x I} Matrix of estimated proportions for each iteration \eqn{i \in \{1, \cdots, I\}} for the first known component}
 \item{pi2}{An \eqn{S_T x I} Matrix of estimated proportions for each iteration \eqn{i \in \{1, \cdots, I\}} for the second known component}
   }
\note{
We develop a three-component deconvolution model, DeMixT, for expression data from a mixture of cancerous tissues, infiltrating immune cells and tumor microenvironment. DeMixT is a frequentist-based method and fast in yielding accurate estimates of cell proportions and compartment-specific expression profiles for two-component and three-component deconvolution problem. Our method promises to provide deeper insight into cancer biomarkers and assist in the development of novel prognostic markers and therapeutic strategies.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\author{
Zeya Wang
}
\references{
}
\keyword{DeMixT.Kernel}
\seealso{
		http://bioinformatics.mdanderson.org/main/DeMix:Overview
 }
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DeMixT}% use one of  RShowDoc("KEYWORDS")
\keyword{DeMixT}% __ONLY ONE__ keyword per line

