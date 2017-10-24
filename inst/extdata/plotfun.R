require(LDheatmap)

gwasvis_helper <- function(gwasdata, block_b, diseaseIdx, pred, 
                           predName="PIP", pred_list, predOnly=FALSE) {
  
  snpid <- gwasdata$snpidByBlocks[[block_b]]
  
  nsnp <- length(snpid)
  
  gwasPval_b <- gwasdata$gwasPval[[block_b]]
  
  causalidx <- gwasdata$is_causal[[block_b]]
  
  if(missing(pred_list)) {
    
    pred_list <- list(gwasPval=gwasPval_b)
    
    if(!missing(pred)) {
      
      if(length(pred) > nsnp) {
        
        pred_list[[predName]] <- pred[snpid,,drop=F]
        
      } else {
        
        pred_list[[predName]] <- pred
      }
    }
  }
  
  D <- ncol(gwasPval_b)
  
  pos <- start(gwasdata$snpGRlist[[block_b]])
  
  df <- do.call(rbind, lapply(1:D, function(d) {
    
    do.call(rbind, lapply(names(pred_list), function(m) {
      
      if(m=="gwasPval") {
        
        data.frame(snpidx=1:nsnp, pos=pos, score=-log10(pred_list[[m]][,d]), 
                   causal=causalidx[,d], disease=sprintf("trait%s",d), method="GWAS -logP")				
      } else {
        data.frame(snpidx=1:nsnp, pos=pos, score=pred_list[[m]][,d], 
                   causal=causalidx[,d], disease=sprintf("trait%s",d), method=m)
      }
    }))
  }))		
  
  if(!missing(diseaseIdx) & D>1) df <- subset(df, disease==diseaseIdx)
  
  gr_sel <- gwasdata$snpGRlist[[block_b]]
  
  chrom <- as.character(seqnames(gr_sel))[1]
  
  locusRange <- sprintf("%s: %s-%s", chrom, min(start(gr_sel)), max(end(gr_sel)))
  
  df$lsnp <- snpid == snpid[which.min(gwasPval_b)]
  
  pos <- start(gr_sel)
  
  names(pos) <- snpid
  
  causal_pos_b <- pos[gwasdata$is_causal[[block_b]]]
  
  if(predOnly) df <- subset(df, method!="GWAS -logP")
  
  gg <- ggplot(df, aes(
    
    # x=snpidx, 
    
    x=pos,
    
    y=score,
    
    color=causal, fill=causal, alpha=lsnp,
    
    shape=lsnp, size=causal
    )) +
    
    geom_point() + theme_bw() +				
    
    scale_color_manual(values=c("blue", "red")) + 
    
    scale_fill_manual(values=c("blue", "red")) + 
    
    scale_alpha_manual(values=c(0.3, 1)) +
    
    scale_shape_manual(values=c(21,23)) +
    
    scale_size_manual(values=c(2, 2)) +
    
    geom_vline(xintercept=causal_pos_b, color="red", alpha=0.2, size=2) +
    
    theme(legend.position="blank") + 
    
    xlab(locusRange) + ylab("Statistical Significance")
  
  # if(missing(diseaseIdx) & (!missing(pred) | !missing(pred_list))) {
  #   gg <- gg + facet_grid(method~disease, scale="free")
  # } 
  
  if(!missing(pred) & !predOnly) {
    gg <- gg + facet_grid(method~., scale="free")
  } 
  
  if(missing(diseaseIdx) & length(pred_list)==1) {
    gg <- gg + ylab("-logP")
  }
  
  if(predOnly) {
    gg <- gg + ylab(predName)
  }
  
  gg
}

gwasvis <- function(gwasdata, block_b, title="", ...) {
  
  pos <- start(gwasdata$snpGRlist[[block_b]])
  
  LDheatmap.addGrob(LDheatmap(gwasdata$hapr2[[block_b]], pos,
                              
                              flip=T, color=heat.colors(20), title=title), 
                    
                    ggplotGrob(gwasvis_helper(gwasdata, block_b, ...)), height=.8)		
}


visld <- function(gwasdata, block_b, title="", ...) {
  
  pos <- start(gwasdata$snpGRlist[[block_b]])
  
  LDheatmap(gwasdata$hapr2[[block_b]], pos, flip=T, color=heat.colors(20))
}

require(ggplot2)

qq <- function(observed, is_flagged) {
  
  expected <- c(1:length(observed)) 
  
  idx <- order(observed)
  observed <- observed[idx]
  
  observed <- -(log10(observed))
  expected <- -(log10(expected / (length(expected)+1)))
  
  if(!missing(is_flagged)) {
    
    qqinp <- data.frame(expected=expected, observed=observed, flagged=is_flagged[idx])
    
    qqp <- ggplot(qqinp, aes(x=expected, y=observed, color=flagged, alpha=flagged)) + 
      geom_point() + stat_binhex() +
      geom_abline(intercept = 0, slope = 1) + theme_bw() +
      scale_color_manual(values=c("blue","red")) + 
      xlab("Expected (-logP)") + ylab("Observed (-logP)") +
      ylim(0, 7) + xlim(0, 7)
    
  } else {
    
    qqinp <- data.frame(expected=expected, observed=observed)
    
    qqp <- ggplot(qqinp, aes(x=expected, y=observed)) + 
      geom_point() + stat_binhex() + 
      geom_abline(intercept = 0, slope = 1) + theme_bw() +
      scale_color_manual(values=c("blue","red")) + 
      xlab("Expected (-logP)") + ylab("Observed (-logP)") +
      ylim(0, 7) + xlim(0, 7)
  }
  
  qqp
}











