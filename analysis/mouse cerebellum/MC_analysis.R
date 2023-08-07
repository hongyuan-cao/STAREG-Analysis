rm(list = ls())
load("./output/MC.RData") # obtained from https://drive.google.com/file/d/13a-GPtPDiwIBNjVbOOOB64x6fQfAh2hn/view?usp=sharing

# Slide-seq
p1 <- sparkx1$res_mtest$combinedPval
names(p1) <- rownames(sparkx2$res_mtest)
# Slide-seqV2
p2 <- sparkx2$res_mtest$combinedPval
names(p2) <- rownames(sparkx2$res_mtest)
alpha = 0.05

##---------------------------------------------------------------------------
## Slide-seq Replicate 2 + Slide-seqV2
##---------------------------------------------------------------------------
overlap <- intersect(names(p1),names(p2))
pvals1 = p1[overlap]
pvals2 = p2[overlap]

library(STAREG)
res.stareg <- stareg(pvals1, pvals2, init.pi0 = FALSE)
genes_rep_stareg <- overlap[which(res.stareg$fdr<=alpha)]

# BH method
p1.bh <- p.adjust(pvals1, method = 'BH')
genes1.bh <- overlap[which(p1.bh<=alpha)]
p2.bh <- p.adjust(pvals2, method = 'BH')
genes2.bh <- overlap[which(p2.bh<=alpha)]
genes_rep_bh <- intersect(genes1.bh, genes2.bh)

# MaxP method
max_pvals <- apply(cbind(pvals1, pvals2), 1, max)
maxp.bh <- p.adjust(max_pvals, method = 'BH')
genes_rep_maxp <- overlap[which(maxp.bh<=alpha)]

# JUMP
library(JUMPrcpp)
jump.obj <- JUMP(pvals1, pvals2, alpha)
jump.thr <- jump.obj$jump.thr
p.max <- jump.obj$p.max
genes_rep_jump <- overlap[which(p.max<=jump.thr)]

## MaRR
source("R/methods.R")
x <- rank(pvals1)
y <- rank(pvals2)
max.rank = apply(cbind(x,y),1,max)
res.marr <- MaRR(max.rank, alpha = alpha)
genes_rep_marr <- overlap[res.marr$which.sig]

## Bogomolov & Heller 2018 (JASA)
library(radjust)
pa = pvals1
pb = pvals2
pa[which(pa==0)] <- 1e-15
pb[which(pb==0)] <- 1e-15
names(pa) <- NULL
names(pb) <- NULL
res.rv18 <- radjust_sym(pa, pb, input_type = "all", directional_rep_claim = FALSE,
                        variant = "adaptive", alpha = alpha)
rv18 <- rep(1, length(pa))
rv18[as.numeric(res.rv18$results_table$name)] <- res.rv18$results_table$r_value
genes_rep_radjust <- overlap[rv18 <= alpha]

save(overlap, pvals1, pvals2, genes_rep_jump, genes_rep_marr, genes_rep_radjust, genes_rep_bh, genes_rep_maxp, genes_rep_stareg, file = "output/MC_results.RData")
