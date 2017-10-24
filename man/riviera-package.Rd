\name{RiVIERA-package}
\alias{RiVIERA-package}
\alias{RiVIERA}
\docType{package}
\title{
  Risk Variant Inference using Epigenomic Reference Annotation to predict Multiple Trait-causing co-localized mutations
}
\description{
  RiVIERA-MT is a full Bayesian framework to fine-map causal variants using GWAS summary statistics in z-scores in one or multiple related traits and large-scale reference annotations in binary or continuous values.
}
\details{
  The goal of RiVIERA is to infer for each SNP in disease their posterior probability of disease association and to detect functional enrichments or depletions from the annotations, taking into account the underlying multi-trait epigenomic covariance. The inputs to the model are:
  
  1. GWAS summary association statistics in one or more traits (Z-scores):
  
  List of matrices, and each matrix is M_d x 1 vector for M_d SNPs in disease d;
  
  2. Linkage disequilibrium (LD) either from the GWAS cohort or reference panel from 1000 Genome consortium:
  
  Nested list of lists. The list contains D items for D traits. Each item is also a list and contain a set of LD matrices for loci of trait d;
  
  3. Functional annotation matrix (binary or continuous) for each SNP:
  
  Same as 1 but contain M_d x K matrix for trait d where K is the number of annotations.
}
\author{
Yue Li

Maintainer: Yue Li <liyue@mit.edu>
}
\references{
  Yue Li and Manolis Kellis. RiVIERA-MT: A Bayesian model to infer risk variants in related traits using summary statistics and functional genomic annotations. bioRxiv 059345; doi: http://dx.doi.org/10.1101/059345
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{RiVIERA}
\examples{
  vignette("RiVIERA")
}
