gwasvis <- function(simdata) {
  
  blockid <- data.frame(
    rsid=unlist(simdata$snpidByBlocks), 
    locusID=unlist(lapply(1:length(simdata$snpidByBlocks), 
                          function(i) rep(i, length(simdata$snpidByBlocks[[i]])))))
  
  inp <- data.frame(rsid=rownames(do.call(rbind, simdata$gwasZval)), 
                    CHR=as.numeric(sub("chr","",seqnames(unlist(simdata$snpGRlist)))),
                    POS=as.numeric(start(unlist(simdata$snpGRlist))), 
                    locusID=blockid$locusID[match(rownames(do.call(rbind,simdata$gwasZval)), blockid$rsid)],
                    Z=do.call(rbind, simdata$gwasZval), 
                    P=do.call(rbind, simdata$gwasPval), 
                    # MAF=simdata$maf, 
                    # N=if(d==1) simdata$sample1 else simdata$sample2,
                    SEGNUMBER=unlist(sapply(1:length(simdata$blocksize), 
                                            function(b) rep(b, simdata$blocksize[b]))),
                    is_causal=unlist(simdata$is_causal))
  
  # plot manhattan plots
  gg <- 
    ggplot(inp, aes(x=POS, y=-log10(P)), color=is_causal) +
    # ggplot(inp, aes(x=POS, y=Z)) +
    geom_point(size=1) + facet_grid(.~locusID, scale="free") + theme_bw() +
    theme(axis.ticks.x=element_blank(), 
          axis.text.x=element_blank(), 
          legend.position="none",
          strip.text.y=element_text(size=10),
          strip.text.y=element_blank(),
          panel.border=element_blank(), 
          strip.background=element_blank(),
          panel.margin.x = unit(0, "lines")) + xlab("") +
    scale_color_manual(values=c("grey", "red"))
  
  gg
}



gwasvis_pipDistrn <- function(simdata, pip_postDistrn, sel) {
  
  locusID <- unlist(lapply(sel, function(l) 
    rep(l, length(simdata$snpidByBlocks[[l]]))))
  
  snpgr <- unlist(simdata$snpGRlist[sel])
  
  pip_postDistrn_df <- as.data.frame(pip_postDistrn)
  
  pip_postDistrn_df$rsid <- unlist(simdata$snpidByBlocks[sel])
  
  pip_postDistrn_df$is_causal <- tru_label
  
  pip_postDistrn_df$locusID <- locusID
  
  pip_postDistrn_df <- melt(pip_postDistrn_df, 
                            id.vars=c("rsid","locusID","is_causal"))
  
  colnames(pip_postDistrn_df)[ncol(pip_postDistrn_df)] <- "pip"
  
  pip_postDistrn_df$variable <- NULL
  
  pip_postDistrn_df$type <- "pip"
  
  pip_postDistrn_df$POS <- as.numeric(start(snpgr))
  
  pip_postDistrn_df$CHR <- as.numeric(sub("chr","",seqnames(snpgr)))
  
  # plot manhattan plots
  gg <- ggplot(pip_postDistrn_df, aes(x=POS, y=-log10(1-pip), color=is_causal)) + 
    geom_boxplot(outlier.size=0.5) + facet_grid(.~locusID, scale="free") + theme_bw() +
    theme(axis.ticks.x=element_blank(), 
          axis.text.x=element_blank(), 
          legend.position="none",
          strip.text.y=element_text(size=10),
          strip.text.y=element_blank(),
          panel.border=element_blank(), 
          strip.background=element_blank(),
          panel.margin.x = unit(0, "lines")) + xlab("") +
    scale_color_manual(values=c("grey", "red"))
  
  gg
}






