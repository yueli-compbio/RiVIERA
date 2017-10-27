---
title: "Risk Variants Inference using Epigenomic Reference Annotations (RiVIERA)"
author: "Yue Li"
date: '`r Sys.Date()`'
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Risk Variants Inference using Epigenomic Reference Annotations (RiVIERA)}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---

## Brief introduction
RiVIERA is a probabilistic framework to infer functional enrichments and prioritize causal variants using summary statistics and epigenomic or other functional genomic annotations. Specifically, RiVIERA can be divided into two stages: (1) RiVIERA-ench: genome-wide enrichment estimations; (2) RiVIERA-fmap: fine-mapping causal variants. We demonstrate how to run both in the following tutorial.

## Installation
You can easily install RiVIERA from Github as follows:

```{r eval=FALSE,echo=FALSE}
library(devtools)

install_github("yueli-compbio/RiVIERA")
```

## Data preparation
To run RiVIERA, we will need the following data. Suppose M SNPs, K annotations, G groups over the K annotations, then the folowing data matrices are expected:

1. GWAS summary statistics in terms p-values in a numerical vector (M x 1)
2. Functional annotation matrix (binary or continuous) for each SNP (M x K)
3. Annotatoin by group binary matrix indicating what group each annotatoin belongs to. This is required only for the group-guided sparse enrichment model (M x K)
4. Linkage disequilibrium (LD) either from the GWAS cohort or reference panel (e.g. 1000 Genome Project). LD matrices are presented in a list with each item representing a square matrix of Pearson correlation between SNPs in that locus. This is required only for the fine-mapping model. 

For illustration purpose, we provide a simulated dataset saved in RData file and loaded as follows. In this simulated dataset, there are 100 causal variants out of 17225 variants and 10 annotations. The causal annotatoins are the first 3 annotatoins.

## Genome-wide enrichment estiamtions
We first demonstrate inferring genome-wide enrichment using an efficient regularized log-likelihood model.

```{r eval=TRUE,echo=TRUE}
library(RiVIERA)

rda <- system.file("extdata/seed10.RData", package="RiVIERA")

if(!file.exists(rda)) {
  rda <- list.files("~/Projects/riviera/simdata",".RData$",full.names=T)[1]  
}

load(rda)

gwas_pval <- unlist(simdata$gwasPval[[1]])

ann <- simdata$ann

res <- riviera_ridge(gwas_pval,
                     ann,
                     ann_w_mu=rep(0, ncol(ann)),
                     minPval=1e-50,
                     alpha0=0.1,
                     pri0=1e-3,
                     epsilon=1,
                     thres=1e-4,
                     max_iter=100)

ann_w <- res$ann_w

rownames(ann_w) <- colnames(ann)

#### significance test ####
scaler <- 1.96 # 95% CI

ann_w_se <- res$ann_w_se[-1] # remove intercept

upper <- ann_w + scaler * ann_w_se
lower <- ann_w - scaler * ann_w_se

ann_w_zscore <- res$ann_w_zscore

enrichment <- data.frame(ann=1:10,
                         causal=c(rep(TRUE,3),rep(FALSE,7)),
                         enrichment_zscore=ann_w_zscore, 
                         ann_w=ann_w, w_se=ann_w_se,
                         ann_w_pval = 2*pnorm(-abs(ann_w_zscore)),
                         upper=upper, lower=lower)

print(enrichment)
```

We now display the enrichment in terms of the estimated p-values, where the causal annotations are disguished from the non-causal ones by the bar color.

```{r, fig.width=6, fig.height=4.8}
library(ggplot2)

ggplot(enrichment, aes(x=ann, y=-log10(ann_w_pval), fill=causal)) +
  geom_bar(stat="identity", position='dodge') + 
  xlab("Annotation") + geom_hline(yintercept = 2) +
  ylab("Model estimates")
```


We can infer posterior probability of associations (PPA) for each SNP being in the risk loci as follows:

```{r eval=TRUE,echo=TRUE}
ppa <- inferPPA(gwas_pval, ann, ann_w, alpha=res$alpha)

rownames(ppa) <- rownames(gwas_pval)

print(head(ppa[order(ppa,decreasing = TRUE),]))
```

The above PPA gives a way infer risk loci. To identify the exact causal variants we will need to account for LD, which we demonstrate next.


## Fine-mapping causal variants
We now demonstrate the fine-mapping component of RiVIERA framework using the same simulated data. Here we incorporate the LD information along with summary statistics z-scores and annotations. Notably, we can also use the above enrichments as the prior information by setting the "ann_w_mu" option. We filter out the loci that don't harbor any causal SNPs and retain 9939 variants and 55 risk loci each harboring at least one causal SNP.

```{r eval=TRUE,echo=TRUE}

risklociIdx <- which(simdata$causalcnt > 0)

riskloci_risd <- unlist(simdata$snpidByBlocks[risklociIdx])

zscore <- simdata$gwasZval[riskloci_risd,,drop=F]

ldmatList <- simdata$transformedLD[risklociIdx]

ann <- simdata$ann[riskloci_risd,]

is_causal <- simdata$is_causal[[1]][riskloci_risd]
```

To fine-map causal variants, we implemented two sampling schemes (configSampleScheme) namely (1) stochastic neighborhood search (nbh); and (2) importance sampling. The former is more accurate when the in-sample LD is used (i.e., LD was estimated directly from the GWAS cohort), which is the case in our simulation. However, the latter has the benifits of being more efficient, emphasizing on the SNPs with large z-scores, and less dependent on the correctness of LD matrix espeically when the LD is estiamted from a reference panel rather than from the GWAS cohort.

```{r eval=TRUE,echo=TRUE}
#### model parameters ####
step_W <- 1e-2
nsteps_W <- 100

step_S <- 1e-3
nsteps_S <- 100

thres <- 1e-5

burninFrac <- 0.2

randomSeed <- 23

set.seed(randomSeed)

# 1: nbh; 2: imp
configSampleScheme <- 1

sampleSize <- 10

max_iter <- 10

res <- riviera_fmap(zscore=zscore,
                    ldmat=ldmatList,
                    ann=ann,
                    ann_w_mu=rep(0,ncol(ann)),
                    ann_w_sigma0=diag(ncol(ann),ncol(ann)),
                    configSampleScheme=configSampleScheme,
                    step_W = step_W,
                    nsteps_W = nsteps_W,
                    step_S = step_S,
                    nsteps_S = nsteps_S,
                    max_iter=max_iter,
                    sampleSize=sampleSize,
                    thres=thres,
                    verboseS=F,
                    verboseW=F,
                    verboseC=F)

burnFrac <- 0.2

burnIdx <- 1:round(burnFrac * res$fit_info$ensembleSize)

# SNP posterior inclusion probabilities (PIP)
pip <- rowMeans(res$ensemble$pip_ensemble[,-burnIdx])

# enrichment weights
w <- res$ensemble$ann_w_ensemble[,-burnIdx,drop=F]

rownames(w) <- sprintf("A%s", 1:nrow(w))

enrichDf <- melt(w)

enrichDf$causal <- enrichDf$X1 %in% sprintf("A%s", 1:3)

annOrder <- levels(enrichDf$X1)[order(as.numeric(sub("A","",levels(enrichDf$X1))))]

enrichDf$X1 <- factor(enrichDf$X1, levels=annOrder)
```

Below we display the distribution of the enrichments and the accuracy of inferring the causal variants in terms of precision-recall.

```{r, fig.width=6, fig.height=4.8}
library(ggplot2)

ggplot(enrichDf, aes(x=X1, y=value, fill=causal)) +
  geom_boxplot() + 
  xlab("Annotation") + 
  ylab("Model estimates")
```


```{r, eval=TRUE,echo=TRUE}
source(system.file("extdata/eval_funs.R", package="RiVIERA"))

library(ROCR)
library(scales)

predList <- list(`RiVIERA`=pip, `-log10(pval)`=-log10(2*pnorm(-abs(zscore))))

powanal_roc <- do.call(rbind, lapply(names(predList), function(x) 
  roceval(predList[[x]], is_causal, x)))

powanal_prc <- do.call(rbind, lapply(names(predList), function(x) 
  prceval(predList[[x]], is_causal, x)))
```

```{r, fig.width=6, fig.height=4.8}
library(gridExtra)

grid.arrange(evalLinePlot(powanal_roc, "ROC"), evalLinePlot(powanal_prc, "PRC"))
```

## Group-guided sparse enrichment learning
The RiVIERA-ridge assumes that the annotations are independent from each other and uses L2-norm (i.e., zeor-mean standard normal prior) to regularize enrichments. However, often the annotations are correlated with each other due to the shared basic regulatory elements. As a result, a naive model will give much interpretable enrichment estimates. Here we address this challenge by harnessing the metadata information about the annotations. 

For instance, the group information can be the group of primariy tissue or cell types. Within each tissue type, there are further refinement on the cell-type-specific annotations. For instance, we can have immune as a general cell group, and under "immune" category, we have different immune cell-types such as T cell, B cells, etc. The basic idea is therefore to iteratively infer the enrichment at the cell-group level and then dissect the cell-type-specific enrichment within the group that exhibit high enrichment.

For illustration purpose, we use Multiple Sclerosis GWAS data as an example. For annotatations, we took 52 non-cell-type specific baseline annotations and 220 cell-type-specific annotatoins that are grouped into 10 cell-groups, namely adrenal/pancreas, cardiovascular, central nervous system (CNS), connective/bone, gastrointestinal, immune, kidney, liver, skeletal muscle, and others (Finucane et al., Nature genetics 2015). These annotations were obtained from the LD score regression database (https://data.broadinstitute.org/alkesgroup/LDSCORE/).

We run our model with L2-norm on the baseline annotatoins and group-lasso prior (L1/L2-norm) on the cell-type-specific annotations. The goal is to find out what cell-type-specific annotations are enriched from Multiple Sclerosis genetic signals. 

Here we impose a sparsity (groupScoreThres) of 0.9, which at each iteration will focus on learning only the annotatoins within the cell group that has enrichment z-score greater than 90% of the maximum z-scores over all of the 10 cell groups.

```{r eval=TRUE,echo=TRUE}

gwasData <- system.file("extdata/Multiple_sclerosis.RData", package="RiVIERA")

if(!file.exists(gwasData)) {
  gwasData <- "~/Projects/haploriv/data/sumstats_genomewide_priorityPrunerOutputs_annot/Multiple_sclerosis.RData"  
}

load(gwasData)

library(openxlsx)

metadata <- system.file("extdata/epiref.xlsx", package="RiVIERA")

if(!file.exists(metadata)) {
  metadata <- "~/Projects/riviera/database/LDSCORE/all/epiref.xlsx"  
}

epiref <- read.xlsx(metadata)

cellgroup <- unique(epiref$cellgroup)

pval <- 2*pnorm(-abs(zscore))

groupMat <- matrix(0, length(cellgroup), nrow(epiref))

for(j in 1:length(cellgroup)) {
  
  groupMat[j,] <- epiref$cellgroup %in% cellgroup[j]
}

groupMat_dummy <- cbind(diag(ncol(ann_bl)+1), matrix(0, ncol(ann_bl)+1, ncol(ann_ct)))

groupMat_ct <- cbind(matrix(0, length(cellgroup), ncol(ann_bl)+1), groupMat)

groupMat <- rbind(groupMat_dummy, groupMat_ct)

annMat <- cbind(1, ann_bl, ann_ct)

# zero-based
targetGroups <- tail(1:nrow(groupMat),length(cellgroup))-1

annGroupNames <- c("dummy", colnames(ann_bl), cellgroup)

targetGroupNames <- annGroupNames[tail(1:nrow(groupMat),length(cellgroup))]

groupScoreThres <- 0.9

res <- riviera_glass(pval=pval,
                     ann=annMat,
                     groupMat=groupMat,
                     targetGroups=targetGroups,
                     targetGroupNames=targetGroupNames,
                     minPval=1e-50,
                     alpha0=0.1,
                     pri0=1e-3,
                     epsilon=0.1,
                     thres=1e-3,
                     max_iter=10,
                     groupScoreThres=groupScoreThres,
                     verbose=1)

K <- ncol(ann_ct)

ann_w <- tail(res$ann_w, K)

rownames(ann_w) <- tail(colnames(annMat), K)

#### significance test ####
scaler <- 1.96 # 95% CI

ann_w_se <- tail(res$ann_w_se[-1],K) # remove intercept

upper <- ann_w + scaler * ann_w_se
lower <- ann_w - scaler * ann_w_se

ann_w_zscore <- tail(res$ann_w_zscore, K)

enrichment <- data.frame(enrichment_zscore=ann_w_zscore, 
                         ann_w=ann_w, w_se=ann_w_se,
                         ann_w_pval = 2*pnorm(-abs(ann_w_zscore)),
                         upper=upper, lower=lower)

rownames(enrichment) <- rownames(ann_w)

enrichW <- cbind(epiref[match(rownames(enrichment), epiref$file_number),], enrichment)
  
enrichW <- enrichW[order(enrichW$enrichment_zscore, decreasing=TRUE),]

print(table(subset(enrichW, ann_w!=0)$cellgroup))

print(head(subset(enrichW, ann_w!=0)))
```

From the above enrichments, we can see that due to the sparsity we impose over the group the annotatoins in the 67 annotatoins in the immune group that the only ones that have non-zero enrichments for MS. This is consistent with the disease biology. After this step, we can also fit the fine-mapping model using the sparse enrichment weights (ann_w) disregarding the weights that are zeros.


## Citation
If you use RiVIERA for your research, please cite the following papers:

Li, Y., Davila-Velderrain, J., & Kellis, M. (2017). A probabilistic framework to dissect functional cell-type-specific regulatory elements and risk loci underlying the genetics of complex traits. bioRxiv. http://doi.org/10.1101/059345








