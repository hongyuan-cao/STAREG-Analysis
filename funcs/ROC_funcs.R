## A revised version of the classical ROC function to make fair comparisons of the power
## We use the empirical FDR, the proportion of false discoveries in all discoveries,
## but not the false positive rate, as x-axis

revisedROC <- function(states, padj){
  fdr = 0
  tpr = 0
  thresholds = c(sort(unique(padj)),1)
  # thresholds = seq(0, 0.2, 0.01)
  
  for (i in 1: length(thresholds)){
    thr = thresholds[i]
    fdr = c(fdr, sum(padj <= thr & !states)/max(sum(padj <= thr), 1))
    tpr = c(tpr, sum(padj <= thr & states)/sum(states))
  }
  # fdr <- round(fdr, 2)
  dups <- duplicated(tpr) & duplicated(fdr)
  fdr <- fdr[!dups]
  tpr <- tpr[!dups]
  
  return(list(fdr = fdr, tpr = tpr))
}

revisedROC2 <- function(states, padj1, padj2){
  fdr = 0
  tpr = 0
  thresholds = c(sort(unique(c(padj1, padj2))),1)
  # thresholds = seq(0, 0.2, 0.01)
  
  for (i in 1: length(thresholds)){
    thr = thresholds[i]
    fdr = c(fdr, sum(padj1 <= thr & padj2 <= thr & !states)/max(sum(padj1 <= thr & padj2 <= thr), 1))
    tpr = c(tpr, sum(padj1 <= thr & padj2 <= thr & states)/sum(states))
  }
  # fdr <- round(fdr, 2)
  dups <- duplicated(tpr) & duplicated(fdr)
  fdr <- fdr[!dups]
  tpr <- tpr[!dups]
  
  return(list(fdr = fdr, tpr = tpr))
}
