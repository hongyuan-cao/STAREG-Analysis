load("./output/MOB.RData") # https://drive.google.com/file/d/1kUr8wW1NheyevfxzP7Fx6Wy3Nndu6j51/view?usp=sharing

##----------------------------------------------------------------------
## Replicability analysis (Replicate 1 + Replicate 8)
##----------------------------------------------------------------------
alpha = 0.05
p1 <- spark1@res_mtest$combined_pvalue
names(p1) <- rownames(spark1@counts)
p2 <- spark8@res_mtest$combined_pvalue
names(p2) <- rownames(spark8@counts)

overlap <- intersect(names(p1),names(p2))
pvals1 = p1[overlap]
pvals2 = p2[overlap]

## STAREG
library(STAREG)
res.repem <- stareg(pvals1, pvals2, init.pi0 = FALSE)
genes_rep_repem <- overlap[which(res.repem$fdr<=alpha)]

## Ad hoc BH method
p1.bh <- p.adjust(pvals1, method = 'BH')
genes1.bh <- overlap[which(p1.bh<=alpha)]
p2.bh <- p.adjust(pvals2, method = 'BH')
genes2.bh <- overlap[which(p2.bh<=alpha)]
genes_rep_bh <- intersect(genes1.bh, genes2.bh)

## MaxP method
max_pvals <- apply(cbind(pvals1, pvals2), 1, max)
maxp.bh <- p.adjust(max_pvals, method = 'BH')
genes_rep_maxp <- overlap[which(maxp.bh<=alpha)]

## MaRR
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

## JUMP
library(JUMP)
jump.obj <- JUMP(pvals1, pvals2, alpha)
jump.thr <- jump.obj$jump.thr
p.max <- jump.obj$p.max
genes_rep_jump <- overlap[p.max <= jump.thr]

save(overlap, pvals1, pvals2, genes_rep_jump, genes_rep_marr, genes_rep_radjust, genes_rep_bh, genes_rep_maxp, genes_rep_repem, file = "output/MOB_results.RData")
