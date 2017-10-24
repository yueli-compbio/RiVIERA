\name{RiVIERA-package}
\alias{RiVIERA-package}
\alias{RiVIERA}
\docType{package}
\title{
  Risk Variant Inference using Epigenomic Reference Annotation
}
\description{
  RiVIERA is a probabilistic framework to infer functional enrichments and prioritize causal variants using summary statistics and epigenomic or other functional genomic annotations. Specifically, RiVIERA can be divided into two stages: (1) RiVIERA-ench: genome-wide enrichment estimations; (2) RiVIERA-fmap: fine-mapping causal variants. We demonstrate how to run both in the following tutorial.
}
\details{
  To run RiVIERA, we will need the following data. Suppose M SNPs, K annotations, G groups over the K annotations, then the folowing data matrices are expected:

1. GWAS summary statistics in terms p-values in a numerical vector (M x 1)

2. Functional annotation matrix (binary or continuous) for each SNP (M x K)

3. Annotatoin by group binary matrix indicating what group each annotatoin belongs to. This is required only for the group-guided sparse enrichment model (M x G)

4. Linkage disequilibrium (LD) either from the GWAS cohort or reference panel (e.g. 1000 Genome Project). LD matrices are presented in a list with each item representing a square matrix of Pearson correlation between SNPs in that locus. This is required only for the fine-mapping model. 
}
\author{
Yue Li

Maintainer: Yue Li <liyue@mit.edu>
}
\references{
  Li, Yue and Davila-Velderrain, Jose and Kellis, Manolis (2017) A probabilistic framework to dissect functional cell-type-specific regulatory elements and risk loci underlying the genetics of complex traits. bioRxiv 10.1101/059345.
}
% Optionally other standard keywords, one per line,
% from the file KEYWORDS in the R documentation.
\keyword{RiVIERA}
\examples{
  vignette("RiVIERA")
}
