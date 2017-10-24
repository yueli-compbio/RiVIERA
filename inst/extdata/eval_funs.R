# ROC
roceval <- function(myscores, labels_true, model) {
  
  pred <- prediction(myscores, labels_true)
  
  perf <- performance( pred, "tpr", "fpr" )
  
  auc <- unlist(slot(performance( pred, "auc" ), "y.values"))
  
  fpr <- unlist(slot(perf, "x.values"))
  
  tpr <- unlist(slot(perf, "y.values"))
  
  cutoffval <- unlist(slot(perf, "alpha.values"))	
  
  rocdf <- data.frame(x= fpr, y=tpr, auc=auc, 
                      
                      cutoff=cutoffval, 
                      
                      evalTypeAUC=sprintf("%s (%s)", model, percent(auc)), 
                      
                      model=model, curveType="ROC")
  
  return(rocdf)
}


# PRC
prceval <- function(myscores, labels_true, model) {
  
  pred <- prediction(myscores, labels_true)
  
  perf <- performance( pred, "prec", "rec" )
  
  rec <- unlist(slot(perf, "x.values"))
  
  prec <- unlist(slot(perf, "y.values"))
  
  cutoffval <- unlist(slot(perf, "alpha.values"))	
  
  prec[is.nan(prec)] <- 0
  
  prec[length(prec)] <- 0
  
  rec[is.nan(rec)] <- 0
  
  auc <- integrate(approxfun(cbind(rec, prec)), lower=0, upper=1,
                   subdivisions=1e5, stop.on.error = FALSE)$value					
  
  prcdf <- data.frame(x=rec, y=prec, auc=auc,
                      evalTypeAUC=sprintf("%s (%s)", model, 
                                          percent(auc)), model=model, curveType="PRC")
  
  return(prcdf)
}


#################### plot ####################
evalLinePlot <- function(mydf, curve, mytitle=NA) {
  
  if(curve=="ROC") {
    x.title <- "False positive rate"	
    y.title <- "True positive rate"
  } else {
    x.title <- "Recall"	
    y.title <- "Precision"		
  }
  
  gg <- ggplot(mydf, aes(x=x, y=y, color=evalTypeAUC)) + 
    
    theme_bw() +
    
    geom_line(aes(linetype=evalTypeAUC), size=0.7) +
    
    scale_y_continuous(y.title, labels=percent) + 		
    
    theme(axis.text.x = element_text(colour = "black"), 
          
          axis.text.y = element_text(colour = "black"),
          
          axis.title.x = element_text(colour = "black"),				
          
          legend.title= element_blank(),	
          
          legend.position=c(0.5,0.3),
          
          plot.title = element_text(colour = "black"),
          
          legend.background = element_blank())
  
  if(!is.na(mytitle)) gg <- gg + ggtitle(mytitle)
  
  if(curve=="ROC") {
    gg + geom_abline(intercept = 0, slope = 1, colour="grey", linetype=2) + 
      scale_x_continuous(x.title, labels=percent)
  } else {
    gg + scale_x_continuous(x.title, labels=percent) # limits=c(0, 0.1)
  }
}	
