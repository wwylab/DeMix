\name{DeMixT.S2}
\alias{DeMixT.S2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Deconvolves expressions of each individual sample for unknown component
}
\description{
This function is designed to estimate the deconvolved expressions of individual mixed tumor samples for unknown component given a subset of genes.
}
\usage{
DeMixT.S2(inputdata, groupid, givenpi,ninteg = 50, filter.option = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{inputdata}{
A matrix of expression data (e.g gene expressions) from reference (e.g. normal) and mixed samples (e.g. mixed tumor samples). It is a \eqn{G x S} matrix where \eqn{G} is the number of genes and \eqn{S} is the number of samples including reference and mixed samples. Samples with the same tissue type should be placed together in columns (e.g. cbind(normal samples, mixed tumor samples)}
  \item{groupid}{
A vector of indicators to denote if the corresponding samples are reference samples or mixed tumor samples. DeMixT is able to deconvolve mixed tumor samples with at most three components. We use 1 and 2 to denote the samples referencing the first and the second known component in mixed tumor samples. We use 3 to indicate mixed tumor samples prepared to be deconvolved. For example, in two-component deconvolution, we have c(1,1,...,3,3) and in three-component deconvolution, we have c(1,1,...,2,2,...,3,3).}
 \item{givenpi}{\eqn{S_{T_N}}-Vector of proportions. Given the number of mixed tumor samples is \eqn{S_T(S_T<S)}, \eqn{S_{T_N}} is set to \eqn{2*S_T} in a three-component setting and \eqn{S_T} in a two-component setting. It is fixed with given proportions for the first and the second known component of mixed tumor samples, or just for one known component when there is just one type of reference tissues. It has the form of Vector \eqn{\pi^1_{N_1}},\eqn{pi^2_{N_1}},...,\eqn{\pi^{S_T}_{N_1}},\eqn{\pi^1_{N_2}},\eqn{\pi^2_{N_2}},...,\eqn{\pi^{S_T}_{N_2}}.}
\item{ninteg}{The number of bins used in numerical integration for computing complete likelihood when expressions are deconvolved. The default value is 50}
\item{filter.option}{The option is used to process or filter zero count from input expression matrix. If it is set to 1, all those genes containing zero count in sample will be removed. If it is set to 2, input data matrix will be added 1 to vanish zero count. This option is used for RNA-seq data. Default is 1.}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{decovExpr}{A matrix of deconvolved expression profiles corresponding to unknown (e.g tumor) component in mixed samples for a given subset of genes. Each row corresponds to one gene and each column corresponds to one sample.}            
 \item{decovMu}{Estimated \eqn{\mu} of log2-normal distribution for unknown component.}
 \item{decovSigma}{Estimated \eqn{\sigma} of log2-normal distribution for unknown component}
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
\keyword{DeMixT.S2}
\seealso{
			http://bioinformatics.mdanderson.org/main/DeMix:Overview
 }
\examples{
#first example given for deconvolving expressions given proportions under the setting of two-component deconvolution
data(simul2)
data(TPi2)
inputdata = as.matrix(simul2)
givenpi <- c(t(as.matrix(TPi2[-2,])))
groupid <-  c(rep(1, 5), rep(3, 10))
res <- DeMixT.S2(inputdata, groupid, givenpi)

#second example given for deconvolving expressions given proportions under the setting of three-component deconvolution
data(simul3)
data(TPi3)
inputdata = as.matrix(simul3)
givenpi <- c(t(as.matrix(TPi3[-3,])))
groupid <-  c(rep(1, 5), rep(2, 5), rep(3, 5))
res <- DeMixT.S2(inputdata, groupid, givenpi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{DeMixT}% use one of  RShowDoc("KEYWORDS")
\keyword{DeMixT}% __ONLY ONE__ keyword per line

