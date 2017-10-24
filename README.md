# RiVIERA
<u>Ri</u>sk <u>V</u>ariant <u>I</u>nference using <u>E</u>pigenomic <u>R</u>eference <u>A</u>nnotation 

## Description
RiVIERA is a probabilistic framework to infer genome-wide enrichments of annotations and leverage the enrichments to fine-map causal variants using GWAS summary statistics and large-scale reference annotations in binary or continuous values. The goal of RiVIERA is to infer for each SNP in disease their posterior probability of disease association and to detect functional enrichments or depletions from the annotations, taking into account the underlying epigenomic covariance.

#### Input
1. GWAS summary association statistics in one or more traits (p-values for enrichments and Z-scores for fine-mapping)
2. Linkage disequilibrium (LD) either from the GWAS cohort or reference panel from 1000 Genome consortium (fine-mapping only)
3. Functional annotation matrix (binary or continuous) for each SNP

#### Important features
1. Efficient genome-wide enrichment inference and group-guided sparse enrichment inference using cell-group level information
2. Bayesian fine-mapping causal SNP and causal annotations allows incorporation of large number of annotations
3. Can model epigenomic covariance of multiple related traits (model complexity is linear to the number of traits)
4. Efficient posteior inference of causal configuration automatically determines the number of causal variants

## Installation
RiVIERA requires Rcpp and RcppArmadillo, which can be easily installed within R via 'install.packages' function.

Download, build and install:

`$ git clone https://github.mit.edu/liyue/RiVIERA`

`$ R CMD build RiVIERA`

`$ R CMD INSTALL RiVIERA_0.9.3.tar.gz`

## Documentation
Function manuals are described in R documents:

`> library(RiVIERA)`

`> ?RiVIERA`

User manual with examples is written in sweave and invoked by vignette:

`> vignette("RiVIERA")`

