# Risk Variants Inference using Epigenomic Reference Annotations (RiVIERA)

## Brief introduction
RiVIERA is a probabilistic framework to infer functional enrichments and prioritize causal variants using summary statistics and epigenomic or other functional genomic annotations. Specifically, RiVIERA can be divided into two stages: (1) RiVIERA-ench: genome-wide enrichment estimations; (2) RiVIERA-fmap: fine-mapping causal variants. We demonstrate how to run both in the following tutorial.

## Genome-wide enrichment estimations


### Data preparation
To run RiVIERA-ench, we will need the following data. Suppose $M$ SNPs, $K$ annotations, $G$ groups over the $K$ annotations, then the folowing data matrices are expected:

1. GWAS summary statistics in terms p-values in a numerical vector ($M\times 1$)
2. Functional annotation matrix (binary or continuous) for each SNP ($M\times K$)
3. Annotatoin by group binary matrix indicating what group each annotatoin belongs to. This is required only for the group-guided sparse enrichment model ($M\times G$)



```r
library(RiVIERA)
```




