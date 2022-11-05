## Generate synthetic data guided by real data
## Compute FDR and Power of all methods across 10 replicates

library(SPARK)
library(RepEM)
library(ggplot2)

source('funcs/SimuFunc.R')

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

methods <- c("MaxP", "BH", "RepEM")

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
  
  # BH
  pvals1.bh <- p.adjust(pvals1, method = "BH")
  pvals2.bh <- p.adjust(pvals2, method = "BH")
  
  # MaxP
  maxp <- apply(cbind(pvals1, pvals2), 1, max)
  pvals.maxp <- p.adjust(maxp, method = "BH")
  
  # RepEM
  res.rep <- RepEM(pvals1, pvals2)
  pvals.rep <- res.rep$fdr.rep
  
  for(j in 1:length(alphas)){
    alpha = alphas[j]
    
    # BH
    res$BH$fdp[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & !truth)/max(sum(pvals1.bh <= alpha & pvals2.bh <= alpha), 1)
    res$BH$pd[i,j] <- sum(pvals1.bh <= alpha & pvals2.bh <= alpha & truth) / sum(truth)
    
    # MaxP
    res$MaxP$fdp[i,j] <- sum(pvals.maxp <= alpha & !truth)/ max(sum(pvals.maxp <= alpha), 1)
    res$MaxP$pd[i,j]  <- sum(pvals.maxp <= alpha & truth) / sum(truth)
    
    # RepEM
    res$RepEM$fdp[i,j] <- sum(pvals.rep <= alpha & !truth)/max(sum(pvals.rep <= alpha), 1)
    res$RepEM$pd[i,j] <- sum(pvals.rep <= alpha & truth) / sum(truth)
  }
}

for (k in 1:length(methods)){
  res[[k]]$fdr <- colMeans(res[[k]]$fdp)
  res[[k]]$pwr <- colMeans(res[[k]]$pd)
}

##------------------------------------------------------------------------------
## FDR plot
##------------------------------------------------------------------------------
fdr.RepEM <- as.data.frame(cbind(alpha = alphas, fdr = res$RepEM$fdr))
fdr.MaxP <- as.data.frame(cbind(alpha = alphas, fdr = res$MaxP$fdr))
fdr.BH <- as.data.frame(cbind(alpha = alphas, fdr = res$BH$fdr))

fdr.BH$method <- "BH"
fdr.MaxP$method <- "MaxP"
fdr.RepEM$method <- "RepEM"

res.fdr <- rbind(fdr.RepEM, fdr.MaxP, fdr.BH)

# 4*4.2 inches
ggplot(res.fdr, aes(alpha, fdr, group = method, color = method, shape = method)) +
  scale_colour_manual(name="",
                      values = c("MaxP"="#6496D2", "BH"="#F4B183", "RepEM"="#8FBC8F")) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0,0.1,0.02)) +
  scale_y_continuous(limits = c(0, 0.16), breaks = seq(0,0.16,0.02)) +
  geom_line(size = 1.5) + xlab("Target FDR level") + ylab("Empirical FDR") +
  theme_bw() +
  geom_abline(intercept = 0, slope = 1, colour = "#44546A", size = 1, linetype = 2) +
  theme(legend.position = "none", # legend.position = c(0.01,1.03)
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.6, linetype = "solid"))

##------------------------------------------------------------------------------
## Power plot
##------------------------------------------------------------------------------
pwr.RepEM <- as.data.frame(cbind(alpha = alphas, pwr = res$RepEM$pwr))
pwr.MaxP <- as.data.frame(cbind(alpha = alphas, pwr = res$MaxP$pwr))
pwr.BH <- as.data.frame(cbind(alpha = alphas, pwr = res$BH$pwr))

pwr.BH$method <- "BH"
pwr.MaxP$method <- "MaxP"
pwr.RepEM$method <- "RepEM"

res.pwr <- rbind(pwr.RepEM, pwr.MaxP, pwr.BH)

library(ggplot2)
# 4*4.2 inches
ggplot(res.pwr, aes(alpha, pwr, group = method, color = method, shape = method)) +
  scale_colour_manual(name="",
                      values = c("MaxP"="#6496D2", "BH"="#F4B183", "RepEM"="#8FBC8F")) +
  scale_x_continuous(limits = c(0, 0.1), breaks = seq(0,0.1,0.02)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0,1,0.1)) +
  geom_line(size = 1.5) + xlab("Target FDR level") + ylab("Power") +
  theme_bw() +
  theme(legend.position = "none",  # legend.position = c(0.01,1.03)
        panel.border = element_blank(),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.x = element_text(size = 15),
        axis.text.y = element_text(size = 15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.6, linetype = "solid"))