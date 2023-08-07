## Generate synthetic data guided by real data
## Compute FDR and Power of all methods across 10 replicates

library(SPARK)
library(STAREG)
library(ggplot2)
library(radjust)
library(JUMP)

source('funcs/SimuFunc.R')
source('funcs/methods.R')

n.rep    = 10
tau1     = 0.2
tau2     = 0.2
ipt      = "Pattern I"
n.gene   = 10000
sig1.str = "Moderate"
sig2.str = "Moderate"
xi00 = 0.85
xi01 = 0.05
xi10 = 0.05
xi11 = 0.05
alphas = seq(0, 0.1, 0.01)

methods <- c("MaxP", "BH", "STAREG", "JUMP", "radjust", "MaRR", "IDR")

data.obj <- list()
res <- list()

for (i in 1: length(methods)){
  res[[methods[i]]]$fdp <- matrix(NA, n.rep, length(alphas))
  res[[methods[i]]]$pd <- matrix(NA, n.rep, length(alphas))
}

for (i in 1: n.rep){
  cat("Replicate ", i, ":\n")
  h = sample(0:3, n.gene, replace = TRUE, prob = c(xi00, xi01, xi10, xi11))
  states1 = rep(0, n.gene)
  states1[which(h==2|h==3)] = 1
  states2 = rep(0, n.gene)
  states2[which(h==1|h==3)] = 1
  data.obj[[i]] <- SimuData(tau1 = tau1,
                            tau2 = tau2,
                            ipt = ipt,
                            sig1.str = sig1.str,
                            sig2.str = sig2.str,
                            truth1   = states1, 
                            truth2   = states2)
  
  pvals1 <- data.obj[[i]]$pvals1
  pvals2 <- data.obj[[i]]$pvals2
  states1 <- data.obj[[i]]$truth1
  states2 <- data.obj[[i]]$truth2
  
  truth = states1 * states2
  
  # Ad hoc BH
  pvals1.bh <- p.adjust(pvals1, method = "BH")
  pvals2.bh <- p.adjust(pvals2, method = "BH")
  
  # MaxP
  maxp <- apply(cbind(pvals1, pvals2), 1, max)
  pvals.maxp <- p.adjust(maxp, method = "BH")

  # MaRR
  rank.p <- rbind(rank(pvals1), rank(pvals2))
  max.rank = apply(cbind(rank.p[1,],rank.p[2,]),1,max)
  
  # STAREG
  res.rep <- stareg(pvals1, pvals2)
  pvals.rep <- res.rep$fdr
  
  for(j in 1:length(alphas)){
    alpha = alphas[j]
    
    # Ad hoc BH
    res$BH$fdp[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & !truth)/max(sum(pvals1.bh <= alpha & pvals2.bh <= alpha), 1)
    res$BH$pd[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & truth) / sum(truth)
    
    # MaxP
    res$MaxP$fdp[i,j] <- sum(pvals.maxp <= alpha & !truth)/ max(sum(pvals.maxp <= alpha), 1)
    res$MaxP$pd[i,j]  <- sum(pvals.maxp <= alpha & truth) / sum(truth)
    
    # STAREG
    res$STAREG$fdp[i,j] <- sum(pvals.rep <= alpha & !truth)/max(sum(pvals.rep <= alpha), 1)
    res$STAREG$pd[i,j] <- sum(pvals.rep <= alpha & truth) / sum(truth)

    # JUMP
    jump.obj <- JUMP(pvals1, pvals2, alpha)
    jump.thr <- jump.obj$jump.thr
    p.max <- jump.obj$p.max
    res$JUMP$fdp[i,j] <- sum(p.max <= jump.thr & !truth)/max(sum(p.max <= jump.thr), 1)
    res$JUMP$pd[i,j] <- sum(p.max <= jump.thr & truth) / sum(truth)
    
    # MaRR
    res.marr <- MaRR(max.rank, alpha = alpha)
    marr.rej <- rep(0, m)
    marr.rej[res.marr$which.sig] = 1
    res$MaRR$fdp[i,j] <- sum(marr.rej & !truth)/ max(sum(marr.rej), 1)
    res$MaRR$pd[i,j] <- sum(marr.rej & truth) / sum(truth)
    
    # Bogomolov & Heller 2018 (JASA)
    pa = pvals1
    pb = pvals2
    pa[which(pa==0)] <- 1e-15
    pb[which(pb==0)] <- 1e-15
    res.rv18 <- radjust_sym(pa, pb, input_type = "all", directional_rep_claim = FALSE,
                            variant = "adaptive", alpha = alpha)
    rv18 <- rep(1, m)
    rv18[as.numeric(res.rv18$results_table$name)] <- res.rv18$results_table$r_value
    res$radjust$fdp[i,j] <- sum(rv18 <= alpha & !truth)/max(sum(rv18 <= alpha), 1)
    res$radjust$pd[i,j] <- sum(rv18 <= alpha & truth) / sum(truth)
  }
}

for (k in 1:length(methods)){
  res[[k]]$fdr <- colMeans(res[[k]]$fdp)
  res[[k]]$pwr <- colMeans(res[[k]]$pd)
}