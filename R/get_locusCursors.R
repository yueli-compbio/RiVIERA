get_locusCursors <- function(hapr2, locus_cnt, zero_based=TRUE) {
  
  if(missing(hapr2) & missing(locus_cnt)) stop("hapr2 or locus_cnt must be provided")
  
  if(!missing(hapr2)) {
    
    locusCursors <- cbind(1:length(hapr2), 1, nrow(hapr2[[1]]))
    
    colnames(locusCursors) <- c("locusID", "SNP_start", "SNP_end")
    
    if(length(hapr2) > 1) {
      for(i in 2:length(hapr2)) {
        locusCursors[i,"SNP_start"] <- locusCursors[i-1,"SNP_start"] + nrow(hapr2[[i-1]])
        locusCursors[i,"SNP_end"] <- locusCursors[i,"SNP_start"] + nrow(hapr2[[i]]) - 1
      }
    }
    
  } else {
    
    locusCursors <- cbind(1:length(locus_cnt), 1, locus_cnt[1])
    
    colnames(locusCursors) <- c("locusID", "SNP_start", "SNP_end")
    
    for(i in 2:length(locus_cnt)) {
      locusCursors[i,"SNP_start"] <- locusCursors[i-1,"SNP_start"] + locus_cnt[i-1]
      locusCursors[i,"SNP_end"] <- locusCursors[i,"SNP_start"] + locus_cnt[i] - 1
    }
  }
  
  if(zero_based) {
    # convert to zero-based indices for C++
    locusCursors[,"SNP_start"] <- locusCursors[,"SNP_start"] - 1
    locusCursors[,"SNP_end"] <- locusCursors[,"SNP_end"] - 1
  }
  
  as.matrix(locusCursors[,-1,drop=F])
}
























