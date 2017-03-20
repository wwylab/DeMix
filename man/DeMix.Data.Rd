\name{DeMix.Data}
\alias{DeMix.Data}
\docType{package}
\title{
  Datasets used in the DeMix R package.
}
\usage{
data(DeMix.Data)
}
\description{
 There are three data objects: simulS1, simulS2, Prop.sim. }

\arguments{
 \item{simulS1}{
A 100 x 20 matrix of simulated RNA-seq count data. It contains 100 genes for expression deconvolution. These first 10 columns correspond to 10 normal samples and second 10 columns correspond to the 10 mixed tumor samples to be deconvolved.
}

 \item{simulS2}{
A 200 x 20 matrix of simulated RNA-seq count data. It contains 200 genes for expression deconvolution. These first 10 columns correspond to 10 normal samples and second 10 columns correspond to the 10 mixed tumor samples to be deconvolved. The dataset simulS1 is a subset of genes in simulS2.
}

\item{Prop.sim}{
A vector, true mixing proportions for samples in both simulS1 and simulS2. }  

\item{simul2}{
	A 120 x 15 matrix of simulated expression data. It contains 120 genes for deconvolution. These first 10 columns correspond to 10 normal samples and second 5 columns correspond to the 5 mixed tumor samples to be deconvolved.
}
\item{simul3}{
	A 120 x 15 matrix of simulated expression data. It contains 120 genes for deconvolution. These first 5 columns correspond to 5 normal samples corresponding to the first known component and those second 5 columns correspond to 5 normal samples corresponding to the second known component, and the last 5 columns correspond to the 5 mixed tumor samples of three component tissues to be deconvolved.
}

\item{cell_line_mix}{
	A 5910 x 20 matrix of normalized RNA-seq count data for mixed cell lines from lung adenocarcinoma, cancer-associated fibroblasts and tumor infiltrating lymphocytes. These first 6 columns correspond to 6 pure samples of cell lines from fibroblasts and those second 6 columns correspond to 6 pure samples of cell lines from lymphocytes, and the last 8 columns correspond to the 8 mixed samples of those three components to be deconvolved.

}
}


\author{
Jaeil Ahn & Zeya Wang
}
\references{
Ahn J, Yuan Y, Parmigiani G, Suraokar MB, Diao L, Wistuba II, and Wang W. Bioinformatics 2013 doi: 10.1093/bioinformatics/btt301.

Wu H, Wang C and Wu Z (2013). “A new shrinkage estimator for dispersion improves differential expression detection in RNA-seq data.” Biostatistics. doi: 10.1093/biostatistics/kxs033.
}
\keyword{DeMix.Data}
\seealso{
 http://bioinformatics.mdanderson.org/main/DeMix:Overview
 }
 
\examples{
 data(simulS1)
}
