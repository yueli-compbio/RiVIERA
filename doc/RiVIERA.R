### R code from vignette source 'RiVIERA.Rnw'

###################################################
### code chunk number 1: loadmypkg
###################################################
library(RiVIERA)


###################################################
### code chunk number 2: simdata
###################################################
load(system.file("extdata/sim_seed103.RData", package="RiVIERA"))
# load("~/Projects/riviera/RiVIERA/inst/extdata/sim_seed103.RData")
sel <- which(simdata$causalcnt > 0)
zscoreList <- list(do.call(rbind, simdata$gwasZval[sel]))
ldmatList <- list(simdata$hapr2[sel])
annList <- list(do.call(rbind, simdata$ann[sel]))
locusCursorList <- list(get_locusCursors(simdata$hapr2[sel]))


###################################################
### code chunk number 3: train
###################################################
set.seed(23)
res <- riviera(
  zscoreList,
  ldmatList,
  annList,
  locusCursorList,
  gwasize=matrix(1e4),
  causalCntPerLocus=1,
  pve_max = 1,
  max_iter =100,
  inferPIP_freq=10,
  samplePVE_freq=10,
  useAnn = TRUE,
  sampleConfig_iter=10,
  step = 0.01,
  nsteps = 100,
  printfreq = 10,
  verbose = T)
burnFrac <- 0.2
burnIdx <- 1:round(burnFrac * length(res$pipList_ensemble))
pred <- averagePIP(res$pipList_ensemble[-burnIdx])[[1]]
# pred <- zscoreList[[1]]
rownames(pred) <- rownames(zscoreList[[1]])


###################################################
### code chunk number 4: vizfm
###################################################
library(ggplot2)
library(ggbio)
library(GenomicRanges)
library(LDheatmap)
library(gridExtra)

source(system.file("extdata/gwasvis.R", package="RiVIERA"))
source(system.file("extdata/plotfun.R", package="RiVIERA"))

ggman_list <- ggann_list <- 
  ggpip_list  <- locusRange_list <- ggbio_tracks_list <- list()

block_list <- c(27,69,23)

for(block_b in block_list) {
  
  ggman <- gwasvis_helper(simdata, block_b)
  
  ggpip <- gwasvis_helper(simdata, block_b, 
                          pred = pred, predOnly = T)
  
  ann_b <- as.data.frame(simdata$ann[[block_b]])
  
  ann_b$rsid <- rownames(ann_b)
  
  pos <- start(simdata$snpGRlist[[block_b]])
  
  names(pos) <- simdata$snpidByBlocks[[block_b]]
  
  causal_pos_b <- pos[names(pos) %in% simdata$causal_snpid]
  
  ann_b$pos <- pos
  
  ann_b <- melt(ann_b, id.vars=c("rsid", "pos"))
  
  ann_b <- subset(ann_b, value > 0)
  
  ann_b$ann[ann_b$value == 1] <- 
    sub(".pval.signal.bigwig.RData","",ann_b$variable[ann_b$value == 1])
  
  gr_sel <- simdata$snpGRlist[[block_b]]
  
  chrom <- as.character(seqnames(gr_sel))[1]
  
  locusRange <- sprintf("%s: %s-%s", chrom, min(start(gr_sel)), max(end(gr_sel)))
  
  ggann <- ggplot(ann_b, aes(x=pos, color=ann, weight=value)) + 
    
    theme_bw() + 
    
    geom_density(stat='bin', position='stack') +
    
    geom_vline(xintercept=causal_pos_b, color="red", alpha=0.2, size=2) +
    
    theme(legend.position="none") +
    
    ylab("Epigenomic activities") + xlab(locusRange)
  
  x <- simdata$hapr2[[block_b]]
  
  dimnames(x) <- list(pos, pos)
  
  ldf <- melt(x)
  
  locusRange_list[[as.character(block_b)]] <- locusRange
  
  ggbio_tracks <- tracks(ggman,
                         ggann,
                         ggpip,
                         xlab=locusRange)
  
  ggbio_tracks_list[[as.character(block_b)]] <- ggbio_tracks
  
  ggman_list[[as.character(block_b)]] <- ggman
  ggann_list[[as.character(block_b)]] <- ggann
  ggpip_list[[as.character(block_b)]] <- ggpip  
}


###################################################
### code chunk number 5: setlocus
###################################################
b <- 1


###################################################
### code chunk number 6: ldm1
###################################################
block_b <- block_list[[b]]
idx <- as.character(block_b)
pos <- start(simdata$snpGRlist[[block_b]])
names(pos) <- simdata$snpidByBlocks[[block_b]]
x <- simdata$hapr2[[block_b]]
dimnames(x) <- list(names(pos),names(pos))
causal_pos_b <- pos[names(pos) %in% simdata$causal_snpid]
LDheatmap(x, pos, flip=T, color=heat.colors(20),
          SNP.name = names(causal_pos_b),
          title="", newpage=TRUE)


###################################################
### code chunk number 7: vizPlots1
###################################################
block_b <- block_list[[b]]
idx <- as.character(block_b)
pos <- start(simdata$snpGRlist[[block_b]])
names(pos) <- simdata$snpidByBlocks[[block_b]]
x <- simdata$hapr2[[block_b]]
dimnames(x) <- list(names(pos),names(pos))
ggbio_tracks_list[[idx]]


###################################################
### code chunk number 8: setlocus
###################################################
b <- 2


###################################################
### code chunk number 9: ldm2
###################################################
block_b <- block_list[[b]]
idx <- as.character(block_b)
pos <- start(simdata$snpGRlist[[block_b]])
names(pos) <- simdata$snpidByBlocks[[block_b]]
x <- simdata$hapr2[[block_b]]
dimnames(x) <- list(names(pos),names(pos))
causal_pos_b <- pos[names(pos) %in% simdata$causal_snpid]
LDheatmap(x, pos, flip=T, color=heat.colors(20),
          SNP.name = names(causal_pos_b),
          title="", newpage=TRUE)


###################################################
### code chunk number 10: vizPlots2
###################################################
block_b <- block_list[[b]]
idx <- as.character(block_b)
pos <- start(simdata$snpGRlist[[block_b]])
names(pos) <- simdata$snpidByBlocks[[block_b]]
x <- simdata$hapr2[[block_b]]
dimnames(x) <- list(names(pos),names(pos))
ggbio_tracks_list[[idx]]


###################################################
### code chunk number 11: setlocus
###################################################
b <- 3


###################################################
### code chunk number 12: ldm3
###################################################
block_b <- block_list[[b]]
idx <- as.character(block_b)
pos <- start(simdata$snpGRlist[[block_b]])
names(pos) <- simdata$snpidByBlocks[[block_b]]
x <- simdata$hapr2[[block_b]]
dimnames(x) <- list(names(pos),names(pos))
causal_pos_b <- pos[names(pos) %in% simdata$causal_snpid]
LDheatmap(x, pos, flip=T, color=heat.colors(20),
          SNP.name = names(causal_pos_b),
          title="", newpage=TRUE)


###################################################
### code chunk number 13: vizPlots3
###################################################
block_b <- block_list[[b]]
idx <- as.character(block_b)
pos <- start(simdata$snpGRlist[[block_b]])
names(pos) <- simdata$snpidByBlocks[[block_b]]
x <- simdata$hapr2[[block_b]]
dimnames(x) <- list(names(pos),names(pos))
ggbio_tracks_list[[idx]]


###################################################
### code chunk number 14: sessi
###################################################
sessionInfo()


