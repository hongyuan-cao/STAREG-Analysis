rm(list = ls())
load("./output/HBC_results.RData")
load("./output/HBC.RData") # obtained from https://drive.google.com/file/d/1QeMJThlvKpqcFBd4Zm9Quuu4VBomdgjV/view?usp=sharing

##----------------------------------------------------------
## Moran's I/Geary's C (spark1: 10X Visium; spark2: ST)
##----------------------------------------------------------
source("./funcs/ValidationFuncs.R")
library(spdep)
stat.res1 <- corrValue(spark1@counts, spark1@location)
stat.res2 <- corrValue(spark2@counts, spark2@location)

genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh)
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap]

MI1 = stat.res1$MI
mi_all <- MI1[overlap]
mi_repem_only <- MI1[genes.repem.only]

factor <- factor(rep(c("STAREG only", "All (10X)"), 
                     times = c(length(genes.repem.only),length(overlap))))
dataset <- data.frame(value = c(mi_repem_only, mi_all), group = factor)
par(mar = c(4, 4.5, 1, 1))
boxplot(value ~ group, dataset, col = c("antiquewhite1", "lightseagreen"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6)

MI2 = stat.res2$MI
mi_all <- MI2[overlap]
mi_repem_only <- MI2[genes.repem.only]

factor <- factor(rep(c("STAREG only", "All (ST)"), 
                     times = c(length(genes.repem.only),length(overlap))))
dataset <- data.frame(value = c(mi_repem_only, mi_all), group = factor)
par(mar = c(4, 4.5, 1, 1))
boxplot(value ~ group, dataset, col = c("antiquewhite1", "lightseagreen"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6)

##--------------------------------------------------------------------
## Summarized spatial expression patterns
##--------------------------------------------------------------------
library(amap)
source("./funcs/PlotFuncs.R")

genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh) 
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap] # SVGs only identified by STAREG 

## Sumaarized patterns for ST data
vst_count <- var_stabilize(spark2@counts) # R function in funcs.R
sig_vst_count <- vst_count[genes_rep_repem, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(spark2@lib_size)))
hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 2
memb <- cutree(hc, k = numC)
# The mean residuals of the three patterns for each location
cent <- NULL
for (k in 1:numC) {
  cent <- cbind(cent, colMeans(sig_vst_res[memb == k, , drop = FALSE]))
}
position_cord <- spark2@location
rownames(position_cord) <- rownames(cent)
# The relative residuals
rel_cent <- t(apply(cent, 1, relative_func))
pd <- setNames(cbind.data.frame(position_cord, rel_cent), c("x", "y", paste0("Pattern ", c("II", "I"))))
MBP <- lapply(1:numC, function(x) {
  pattern_plot2(pd, x, xy = T, main = T, titlesize = 1.6)
})
grid.arrange(grobs = MBP[numC:1], nrow = 2) # 6 * 3 inches

## Summarized patterns for 10X Visium data
vst_count <- var_stabilize(spark1@counts) # R function in funcs.R
sig_vst_count <- vst_count[genes_rep_repem, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(spark1@lib_size)))

hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 3
memb <- cutree(hc, k = numC)
# The mean residuals of the three patterns for each location
cent <- NULL
for (k in 1:numC) {
  cent <- cbind(cent, colMeans(sig_vst_res[memb == k, , drop = FALSE]))
}
position_cord <- spark1@location
rownames(position_cord) <- rownames(cent)
# The relative residuals
rel_cent <- t(apply(cent, 1, relative_func))
pd <- setNames(cbind.data.frame(position_cord, rel_cent), c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP <- lapply(1:numC, function(x) {
  pattern_plot2(pd, x, xy = T, main = T, titlesize = 3.5)
})
grid.arrange(grobs = MBP[numC:1], nrow = 2) # 11.5 * 11.5 inches

##--------------------------------------------------------------------
## UMAP clusters
##--------------------------------------------------------------------
library(umap)
library(amap)
source("./funcs/PlotFuncs.R")

## ST result
vst_count <- var_stabilize(spark2@counts) # R function in funcs.R
sig_vst_count <- vst_count[genes_rep_repem, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(spark2@lib_size)))

hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 2
memb <- cutree(hc, k = numC)

counts <- spark2@counts[genes_rep_repem,]
counts <- as.matrix(counts)
scaled_counts <- t(scale(t(counts)))
umap_results <- umap::umap(scaled_counts)

labels <- memb
labels[which(labels == 1)] <- "I"
labels[which(labels == 2)] <- "II"

df <- data.frame(umap_results$layout[,1], umap_results$layout[,2], labels)
colnames(df) <- c("UMAP1","UMAP2", "Cluster")
ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster)) + geom_point() +
  scale_colour_manual(name="Cluster",  
                      values = c("I"="#F4A460", "II"="Steelblue"))


## 10X Visium result
vst_count <- var_stabilize(spark1@counts) # R function in funcs.R
sig_vst_count <- vst_count[genes_rep_repem, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(spark1@lib_size)))
hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 2
memb <- cutree(hc, k = numC)

counts <- spark1@counts[genes_rep_repem,]
counts <- as.matrix(counts)
scaled_counts <- t(scale(t(counts)))
umap_results <- umap::umap(scaled_counts)

labels <- memb
labels[which(labels == 1)] <- "I"
labels[which(labels == 2)] <- "II"

df <- data.frame(umap_results$layout[,1], umap_results$layout[,2], labels)
colnames(df) <- c("UMAP1","UMAP2", "Cluster")
ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster)) + geom_point() +
  scale_colour_manual(name="Cluster",  
                      values = c("I"="#F4A460", "II"="Steelblue"))

##--------------------------------------------------------------------
## Spatial expression patterns of randomly selected genes
##--------------------------------------------------------------------
source("./funcs/PlotFuncs.R")
genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh) 
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap]

randi <- sample(length(genes.repem.only), 20, replace = FALSE)
gene_plot <- genes.repem.only[randi]

## ST
vst_ct <- var_stabilize(spark2@counts) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)
pltdat <- cbind.data.frame(spark2@location[,1:2],rel_vst_ct)
genetitle <- gene_plot
pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot2(pltdat,x,main=T,titlesize=2.2,title=genetitle[x])})
grid.arrange(grobs=pp, nrow=2)


## 10X Visium
vst_ct <- var_stabilize(spark1@counts) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)
pltdat <- cbind.data.frame(spark1@location[,1:2],rel_vst_ct)
genetitle <- gene_plot
pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot2(pltdat,x,main=T,titlesize=4.8,title=genetitle[x])})
grid.arrange(grobs=pp, nrow=2)

##--------------------------------------------------------------------
## Validation gene sets
##--------------------------------------------------------------------
# Genes related to breast cancer in Harmonizome database
library(rjson)
bh.only    <- genes_rep_bh[!genes_rep_bh%in%genes_rep_maxp]
repem.only <- genes_rep_repem[!genes_rep_repem%in%genes_rep_maxp]

Harmonizome_DISEASES <- fromJSON(paste(readLines("./validation/human breast cancer/Harmonizome_DISEASES.json")))
Harmonizome_GAD      <- fromJSON(paste(readLines("./validation/human breast cancer/Harmonizome_GAD.json")))
Harmonizome_GWAS     <- fromJSON(paste(readLines("./validation/human breast cancer/Harmonizome_GWAS.json")))
Harmonizome_OMIM     <- fromJSON(paste(readLines("./validation/human breast cancer/Harmonizome_OMIM.json")))
Harmonizome_PSP      <- fromJSON(paste(readLines("./validation/human breast cancer/Harmonizome_PSP.json")))
Harmonizome_GWASdb   <- fromJSON(paste(readLines("./validation/human breast cancer/Harmonizome_GWASdb.json")))

Harmonizome_DISEASES <- Harmonizome_DISEASES[["associations"]]
Harmonizome_GAD      <- Harmonizome_GAD[["associations"]]
Harmonizome_GWAS     <- Harmonizome_GWAS[["associations"]]
Harmonizome_OMIM     <- Harmonizome_OMIM[["associations"]]
Harmonizome_PSP      <- Harmonizome_PSP[["associations"]]
Harmonizome_GWASdb     <- Harmonizome_GWASdb[["associations"]]

Harmonizome.genes <- NULL
for (i in 1:length(Harmonizome_DISEASES)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_DISEASES[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_GAD)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_GAD[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_GWAS)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_GWAS[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_OMIM)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_OMIM[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_PSP)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_PSP[[i]][["gene"]][["symbol"]])
}
for (i in 1:length(Harmonizome_GWASdb)){
  Harmonizome.genes   <- c(Harmonizome.genes, Harmonizome_GWASdb[[i]][["gene"]][["symbol"]])
}

Harmonizome.genes <- Harmonizome.genes[!duplicated(Harmonizome.genes)] # 3505 genes

length(intersect(genes_rep_maxp, Harmonizome.genes)) # 216/500
length(intersect(genes_rep_bh, Harmonizome.genes)) # 217/504
length(intersect(genes_rep_repem, Harmonizome.genes)) # 241/566

length(intersect(bh.only, Harmonizome.genes)) # 1/4
length(intersect(repem.only, Harmonizome.genes)) # 25/66

# Differential expressed genes of subsets within major cell lineages (Wu et al., 2021)
library(readxl)

sc.genes <- read_excel("./validation/human breast cancer/Differentially expressed genes of subsets within major cell lineages.xlsx")
sc.genes <- unique(sc.genes$gene)

length(intersect(genes_rep_maxp, sc.genes)) # 372/500
length(intersect(genes_rep_bh, sc.genes)) # 374/504
length(intersect(genes_rep_repem, sc.genes)) # 418/566

length(intersect(bh.only, sc.genes)) # 2/4
length(intersect(repem.only, sc.genes)) # 46/66

##--------------------------------------------------------------------
## GO enrichment analysis
##--------------------------------------------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(ggrepel)
library(dplyr)
# options(connectionObserver = NULL) # run if library(org.Mm.eg.db) failed

genes.all <- bitr(overlap, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
## GO enrichment of BH
genes.bh <- bitr(genes_rep_bh, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
go.bh <- enrichGO(gene          = genes.bh$ENTREZID,
                  universe      = genes.all$ENTREZID,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "All" ,
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.99,
                  qvalueCutoff  = 0.99,
                  readable      = TRUE,
                  pool=TRUE)
sum(go.bh$p.adjust<0.05)  # 532

## GO enrichment of STAREG
genes.repem <- bitr(genes_rep_repem, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
go.repem <- enrichGO(gene          = genes.repem$ENTREZID,
                     universe      = genes.all$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "All" ,
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.99,
                     qvalueCutoff  = 0.99,
                     readable      = TRUE,
                     pool=TRUE)
sum(go.repem$p.adjust<0.05) # 554
go.overlap   <- intersect(go.bh@result[go.bh@result$p.adjust<0.05,]$ID, go.repem@result[go.repem@result$p.adjust<0.05,]$ID)
go.repem.sig <- go.repem[go.repem@result$p.adjust<0.05,]
go.repem.only <- go.repem.sig[!go.repem.sig$ID%in%go.overlap,]

write.csv(go.repem,file = "./output/HBC_go_repem.csv",quote = FALSE)
write.csv(go.repem.only,file = "./output/HBC_go_repem_only.csv",quote = FALSE)

results <- go.repem@result
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
  dplyr::select(-len) %>%
  
  # Add this info to the initial dataset
  left_join(results, ., by=c("ONTOLOGY"="ONTOLOGY")) %>%
  
  # Add a cumulative position of group
  arrange(ONTOLOGY, BP) %>%
  mutate(BPcum=BP+tot)

axisdf = don %>% group_by(ONTOLOGY) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
annotated = c("cell adhesion", "focal adhesion", "peptidase inhibitor activity",
              "extracellular matrix", "extracellular matrix structural constituent",
              "collagen metabolic process", "B cell mediated immunity", "regulation of leukocyte activation",
              "lymphocyte mediated immunity")

ggplot(don, aes(x = BPcum, y=-log10(pvalue))) +
  geom_point(aes(color = as.factor(ONTOLOGY), size = Count), alpha=0.8) +
  scale_colour_manual(name="",
                      values = c("BP"="Tan", "CC"="DarkSeaGreen", "MF"="#6496D2")) +
  scale_x_continuous(label = axisdf$ONTOLOGY, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 35)) +     # remove space between plot area and x axis
  geom_hline(yintercept = -log10(0.0089), color = '#545454', size= 1.2, linetype = "dashed") +
  guides(color = "none") + theme_bw() +
  theme(
    legend.position = c(0.02,0.98),
    legend.justification = c(0.02,0.98),
    panel.border = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 12), # face = "bold",
    axis.text.y = element_text(size = 12),
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

genes.all <- bitr(overlap, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
genes.bh <- bitr(genes_rep_bh, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)

kegg.bh <- enrichKEGG(gene         = genes.bh$ENTREZID,
                      organism     = 'hsa',
                      universe     = genes.all$ENTREZID,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
sum(kegg.bh$p.adjust<0.05) # 21

genes.repem <- bitr(genes_rep_repem, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Hs.eg.db)
kegg.repem <- enrichKEGG(gene         = genes.repem$ENTREZID,
                      organism     = 'hsa',
                      universe     = genes.all$ENTREZID,
                      pvalueCutoff = 0.9,
                      qvalueCutoff =0.9)
sum(kegg.repem$p.adjust<0.05) # 25

sig.bh <- filter(kegg.bh, p.adjust<.05)
sig.repem <- filter(kegg.repem, p.adjust<.05)
kegg.overlap    <- intersect(sig.bh@result$ID, sig.repem@result$ID)
kegg.repem.only  <- sig.repem[!sig.repem@result$ID%in%kegg.overlap]

sig.repem.only <- sig.repem %>% filter(ID %in% kegg.repem.only$ID)
kegg.bar <- barplot(sig.repem,color = "pvalue")
kegg.dot <- dotplot(sig.repem.only,color = "pvalue") + 
  theme( 
    legend.position = c(0.98,0.02),
    legend.justification = c(0.98,0.02),
    panel.grid.major = element_blank(),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 12), # face = "bold",
    axis.title = element_text(size = 15),
    panel.grid.minor = element_blank() 
  ) + scale_color_continuous()
