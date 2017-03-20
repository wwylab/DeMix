# DeMix

Deconvolution models for mixed transcriptomes from heterogeneous tumor samples with two or three components using expression data from RNAseq or microarray platforms.

In the DeMixBayes, the deconvolution is designed for a probeset. We focus on the estimation of pi, the purity for the unknown tumor component. The following three steps comprise the selection procedure of this subset.
Step 1: The probes are selected with small standard deviations of the log2- transformed intensity for normal tissues from historical samples. Step 2: Those selected probes are sorted by their standard deviations of raw intensity and those
ones with the largest standard deviations are prepared for the next step. Step 3: The DeMix is run on all probes to estimate mu and pi for each probeset.

In the DeMixT, it is designed to finish the whole pipeline of deconvolution in a setting of two or three components. DeMixT.S1 function is designed to estimate the proportions of all mixed samples for each mixture component with or without subsetting gene set. DeMixT.S2 function is designed to estimate the deconvolved expressions of individual mixed tumor samples for unknown component given a subset of genes.
