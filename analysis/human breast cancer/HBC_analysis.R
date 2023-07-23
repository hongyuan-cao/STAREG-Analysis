load("./output/HBC.RData") # obtained from https://drive.google.com/file/d/1QeMJThlvKpqcFBd4Zm9Quuu4VBomdgjV/view?usp=sharing

##-----------------------------------------------------------------------
## Replicability analysis
##-----------------------------------------------------------------------
alpha = 0.05
p1 <- spark1@res_mtest$combined_pvalue
names(p1) <- rownames(spark1@counts)
p2 <- spark2@res_mtest$combined_pvalue
names(p2) <- rownames(spark2@counts)

overlap <- intersect(names(p1),names(p2))
pvals1 = p1[overlap]
pvals2 = p2[overlap]


library(STAREG)
rep.obj <- stareg(pvals1, pvals2)
genes_rep_repem <- overlap[which(rep.obj$fdr<=alpha)]

# BY method
p1.by <- p.adjust(pvals1, method = 'BY')
genes1.by <- names(pvals1)[which(p1.by<=alpha)]
p2.by <- p.adjust(pvals2, method = 'BY')
genes2.by <- names(pvals2)[which(p2.by<=alpha)]
genes_rep_by <- intersect(genes1.by, genes2.by)

# BH method
p1.bh <- p.adjust(pvals1, method = 'BH')
genes1.bh <- names(pvals1)[which(p1.bh<=alpha)]
p2.bh <- p.adjust(pvals2, method = 'BH')
genes2.bh <- names(pvals2)[which(p2.bh<=alpha)]
genes_rep_bh <- intersect(genes1.bh, genes2.bh)

# MaxP method
max_pvals <- apply(cbind(pvals1, pvals2), 1, max)
maxp.bh <- p.adjust(max_pvals, method = 'BH')
genes_rep_maxp <- overlap[which(maxp.bh<=alpha)]

# Venn diagram
library(VennDiagram)
venn.diagram(
  x = list(genes_rep_repem, genes_rep_bh, genes_rep_maxp),
  category.names = c("STAREG", "BH", "MaxP"),
  fill = c("#E69F00", "coral", "#63769d"),
  filename = 'output/venn_HBC.tiff',
  output = TRUE,
  margin = 0.05,
  cex = 1.5,
  cat.cex = 1.5
)

save(overlap, pvals1, pvals2, genes_rep_bh, genes_rep_maxp, genes_rep_repem, file = "output/HBC_results.RData")
