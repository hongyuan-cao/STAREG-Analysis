rm(list = ls())
load("./output/MC.RData") # obtained from https://drive.google.com/file/d/13a-GPtPDiwIBNjVbOOOB64x6fQfAh2hn/view?usp=sharing
load("./output/MC_results.RData")

##----------------------------------------------------------
## Moran's I
##----------------------------------------------------------
source("./funcs/ValidationFuncs.R")

## SVGs identified by STAREG but not by BH
genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh) # overlap between STAREG and BH
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap]

# Slide-seq
stat.res1 <- corrValue(counts1, location1[,2:3])
MI1 <- stat.res1$MI
mi_all1 <- MI1[overlap]
mi_repem_only1 <- MI1[genes.repem.only]

factor <- factor(rep(c("STAREG only", "All (Slide-seq)"),
                     times = c(length(genes.repem.only),length(overlap))))
dataset1 <- data.frame(value = c(mi_repem_only1, mi_all1), group = factor)
par(mar = c(4, 4.5, 1, 1))
boxplot(value ~ group, dataset1, col = c("antiquewhite1", "lightseagreen"),
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6)

# Slide-seqV2
stat.res2 <- corrValue(counts2, location2)
MI2 <- stat.res2$MI
mi_all2 <- MI2[overlap]
mi_repem_only2 <- MI2[genes.repem.only]

factor <- factor(rep(c("STAREG only", "All (Slide-seqV2)"),
                     times = c(length(genes.repem.only),length(overlap))))
dataset2 <- data.frame(value = c(mi_repem_only2, mi_all2), group = factor)
par(mar = c(4, 4.5, 1, 1))
boxplot(value ~ group, dataset2, col = c("antiquewhite1", "lightseagreen"),
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6)

##--------------------------------------------------------------------
## Summarized spatial expression patterns
##--------------------------------------------------------------------
library(amap)
source("./funcs/PlotFuncs.R")

# Slide-seq
vst_ct1 <- var_stabilize(counts1[overlap,]) 
sig_vst_count1 <- vst_ct1[genes_rep_repem, ]
sig_vst_res1 <- t(apply(sig_vst_count1, 1, LMReg, T = log(apply(counts1[overlap,], 2, sum))))
hc1 <- hcluster(sig_vst_res1, method = "euc", link = "ward", nbproc = 1,
                doubleprecision = TRUE)
numC <- 3
memb1 <- cutree(hc1, k = numC)
# The mean residuals of the three patterns for each location
cent1 <- NULL
for (k in 1:numC) {
  cent <- cbind(cent1, colMeans(sig_vst_res1[memb1 == k, , drop = FALSE]))
}
location1_rev <- location1
location1_rev$y <- -location1$y
position_cord1 <- location1_rev[,2:3]
rownames(position_cord1) <- rownames(cent1)
# The relative residuals
rel_cent1 <- t(apply(cent1, 1, relative_func))
pd <- setNames(cbind.data.frame(position_cord1, rel_cent1), c("x", "y", paste0("Pattern ", c("I", "II", "III"))))
MBP1 <- lapply(1:numC, function(x) {
  pattern_plot3(pd, x, xy = T, main = T, titlesize = 12.8)
})
grid.arrange(grobs = MBP1[1:3], nrow = 3)


## Slide-seqV2
vst_ct2 <- var_stabilize(counts2[overlap,]) # R function in funcs.R
sig_vst_count2 <- vst_ct2[genes_rep_repem, ]
sig_vst_res2 <- t(apply(sig_vst_count2, 1, LMReg, T = log(apply(counts2[overlap,], 2, sum))))

hc2 <- hcluster(sig_vst_res2, method = "euc", link = "ward", nbproc = 1,
                doubleprecision = TRUE)
numC <- 3
memb2 <- cutree(hc2, k = numC)
# The mean residuals of the three patterns for each location
cent2 <- NULL
for (k in 1:numC) {
  cent2 <- cbind(cent2, colMeans(sig_vst_res2[memb2 == k, , drop = FALSE]))
}
location2_rev <- location2
location2_rev$y <- -location2_rev$y
position_cord2 <- location2_rev
rownames(position_cord2) <- rownames(cent2)
# The relative residuals
rel_cent2 <- t(apply(cent2, 1, relative_func))
rel_cent2 <- rel_cent2[,c(2,1,3)]
pd2 <- setNames(cbind.data.frame(position_cord2, rel_cent2), c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP2 <- lapply(1:numC, function(x) {
  pattern_plot3(pd2, x, xy = T, main = T, titlesize = 9.5)
})
grid.arrange(grobs = MBP2[numC:1], nrow = 3)

##--------------------------------------------------------------------
## Spatial expression patterns of randomly selected genes
##--------------------------------------------------------------------
source("./funcs/PlotFuncs.R")

randi <- sample(length(genes.repem.only), 10, replace = FALSE)
gene_plot <- genes.repem.only[randi]
# "Ppp1r17", "Luc7l3", "Atp6v1g2", "Chchd2", "2010107E04Rik", "App", "Malat1", "Cox8a", "Mtss1", "Meg3"
# gene_plot <- c("Meg3", "Ppp1r17") # two representative genes for the main patterns

## Slide-seq
vst_ct1 <- var_stabilize(counts1[overlap,]) 
sig_vst_ct1 <- vst_ct1[gene_plot, ]
rel_vst_ct1 <- apply(sig_vst_ct1, 1, relative_func)
location1_rev <- location1
location1_rev$y <- -location1$y
pltdat <- cbind.data.frame(location1_rev[,2:3],rel_vst_ct1)
genetitle <- gene_plot
pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot3(pltdat,x,main=T,titlesize=12.8,title=genetitle[x])})
grid.arrange(grobs=pp, nrow=2)

## Slide-seqV2
vst_ct2 <- var_stabilize(counts1[overlap,]) # R function in funcs.R
sig_vst_ct <- vst_ct2[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)
location2_rev <- location2
location2_rev$y <- -location2_rev$y
pltdat <- cbind.data.frame(location2_rev,rel_vst_ct)
genetitle <- gene_plot
pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot3(pltdat,x,main=T,titlesize=10.5,title=genetitle[x])})
grid.arrange(grobs=pp, nrow=2)

##--------------------------------------------------------------------
## Mouse cerebellum related studies validation
##--------------------------------------------------------------------
library(readxl)
bh.only    <- genes_rep_bh[!genes_rep_bh%in%genes_rep_maxp]
repem.only <- genes_rep_repem[!genes_rep_repem%in%genes_rep_maxp]

## Genes related to the cerebellum in Harmonizome database (2431)
library(rjson)
library(stringr)

# BioGPS mouse cell type
Harmonizome_BioGPS <- fromJSON(paste(readLines("./validation/mouse cerebellum/Harmonizome_BioGPS.json")))
Harmonizome_BioGPS <- Harmonizome_BioGPS[["associations"]]
Harmonizome.genes <- NULL

for (i in 1:length(Harmonizome_BioGPS)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_BioGPS[[i]][["gene"]][["symbol"]])
}

Harmonizome.genes <- Harmonizome.genes[!duplicated(Harmonizome.genes)]
Harmonizome.genes <- tolower(Harmonizome.genes)
Harmonizome.genes <- str_to_title(Harmonizome.genes) 

length(intersect(genes_rep_maxp, Harmonizome.genes)) # 89/286
length(intersect(genes_rep_bh, Harmonizome.genes)) # 137/476
length(intersect(genes_rep_repem, Harmonizome.genes)) # 231/844

length(intersect(bh.only, Harmonizome.genes)) # 48/190
length(intersect(repem.only, Harmonizome.genes)) # 142/558


## Kozareva et al.
cluster_genes <- read_xlsx("./validation/mouse cerebellum/Kozareva et al.xlsx", col_names = TRUE)
cluster_genes = cluster_genes[which(abs(cluster_genes$logFC)>=1.5),] 
cluster_genes = as.vector(unique(cluster_genes$gene)) # 955

length(intersect(genes_rep_maxp, cluster_genes)) # 85/286
length(intersect(genes_rep_bh, cluster_genes)) # 117/476
length(intersect(genes_rep_repem, cluster_genes)) # 185/844

length(intersect(bh.only, cluster_genes)) # 32/190
length(intersect(repem.only, cluster_genes)) # 100/558

##--------------------------------------------------------------------
## GO enrichment
##--------------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
# options(connectionObserver = NULL) # run if library(org.Mm.eg.db) failed

genes.repem.only <- genes_rep_repem[!genes_rep_repem%in%genes_rep_bh]

genes.all <- bitr(overlap, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh) # overlap between STAREG and BH
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap]
genes_repem_only <- bitr(genes.repem.only, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)

go.repem.only <- enrichGO(gene          = genes_repem_only$ENTREZID,
                          universe      = genes.all$ENTREZID,
                          OrgDb         = org.Mm.eg.db,
                          ont           = "All" ,
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.99,
                          qvalueCutoff  = 0.99,
                          readable      = TRUE,
                          pool=TRUE)
sum(go.repem.only$p.adjust<0.05) # 271
sig.go <- filter(go.repem.only, p.adjust<.05)
write.csv(go.repem.only,file = "./output/MC_go_repem_only.csv",quote = FALSE)

results <- go.repem.only@result
results <- results[sample(nrow(results)),]
results <- results[,c("ID","ONTOLOGY", "pvalue", "Count", "Description")]
results[,"Category"] = NA
results$Category[which(results$ONTOLOGY == "BP")] = 1
results$Category[which(results$ONTOLOGY == "CC")] = 2
results$Category[which(results$ONTOLOGY == "MF")] = 3
results <- results[order(results$Category),]

results[,"BP"] <- NA
results$BP[which(results$Category == 1)] = 1:sum(results$Category == 1)
results$BP[which(results$Category == 2)] = 1:sum(results$Category == 2)
results$BP[which(results$Category == 3)] = 1:sum(results$Category == 3)

don <- results %>%
  
  # Compute group size
  group_by(ONTOLOGY) %>%
  summarise(len=max(BP)) %>%
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(len)-len) %>%
  select(-len) %>%
  
  # Add this info to the initial dataset
  left_join(results, ., by=c("ONTOLOGY"="ONTOLOGY")) %>%
  
  # Add a cumulative position of group
  arrange(ONTOLOGY, BP) %>%
  mutate(BPcum=BP+tot)

axisdf = don %>% group_by(ONTOLOGY) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
annotated = c("neuron to neuron synapse", "modulation of chemical synaptic transmission", "regulation of trans-synaptic signalling",
              "regulation of membrane potential", "synaptic vesicle", "neurotransmitter transport",
              "signal release from synapse", "regulation of synaptic plasticity", "structural constituent of synapse")


ggplot(don, aes(x = BPcum, y=-log10(pvalue))) +
  geom_point(aes(color = as.factor(ONTOLOGY), size = Count), alpha=0.8) +
  scale_colour_manual(name="",
                      values = c("BP"="Tan", "CC"="DarkSeaGreen", "MF"="#6496D2")) +
  scale_x_continuous(label = axisdf$ONTOLOGY, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 13.5)) +     # remove space between plot area and x axis
  geom_hline(yintercept = -log10(0.0037), color = '#545454', size= 1.2, linetype = "dashed") +
  guides(color = "none") +
  theme_bw() +
  theme(
    legend.position = c(0.02,0.98),
    legend.justification = c(0.02,0.98),
    panel.border = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(size = 12.5),
    axis.text.y = element_text(size = 11),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, linetype = "solid")
  ) +
  geom_text_repel(
    data = don[don$Description %in% annotated,],
    aes(label = Description),
    size = 4,
    segment.color = "black", show.legend = FALSE)

##--------------------------------------------------------------------
## KEGG enrichment analysis
##--------------------------------------------------------------------
# install.packages('R.utils')
R.utils::setOption( "clusterProfiler.download.method",'auto')
library(cowplot)

genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh) # overlap between STAREG and BH
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap]
genes_repem_only <- bitr(genes.repem.only, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
kegg.repem.only <- enrichKEGG(gene         = genes_repem_only$ENTREZID,
                              organism     = 'mmu',
                              universe     = genes.all$ENTREZID,
                              pvalueCutoff = 0.9,
                              qvalueCutoff =0.9)
sum(kegg.repem.only$p.adjust<0.05) # 24

sig.kegg <- filter(kegg.repem.only, p.adjust<.05)
dim(sig.kegg)
kegg.bar <- barplot(sig.kegg,showCategory=20,color = "pvalue")
kegg.dot <- dotplot(sig.kegg,showCategory=20,color = "pvalue") + 
  theme( 
    legend.position = c(0.98,0.02),
    legend.justification = c(0.98,0.02),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() 
  ) + scale_color_continuous()



