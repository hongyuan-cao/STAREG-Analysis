rm(list = ls())
load("./output/MOB_results.RData") 
load("./output/MOB.RData") # obtained from https://drive.google.com/file/d/1kUr8wW1NheyevfxzP7Fx6Wy3Nndu6j51/view?usp=sharing

##----------------------------------------------------------
## Moran's I/Geary's C: count in Rep1 and Rep8
##----------------------------------------------------------
source("./funcs/ValidationFuncs.R")
library(spdep)

genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh) 
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap]

# Rep1
stat.res1 <- corrValue(counts1, location1)
MI1 <- stat.res1$MI
mi_all1 <- MI1[overlap]
mi_repem_only1 <- MI1[genes.repem.only]

factor <- factor(rep(c("STAREG only", "All (Rep1)"), 
                     times = c(length(genes.repem.only),length(overlap))))
dataset1 <- data.frame(value = c(mi_repem_only1, mi_all1), group = factor)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset1, col = c("antiquewhite1", "lightseagreen"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6)

# Rep8
stat.res2 <- corrValue(counts8, location8)
MI2 <- stat.res2$MI
mi_all2 <- MI2[overlap]
mi_repem_only2 <- MI2[genes.repem.only]

factor <- factor(rep(c("STAREG only", "All (Rep8)"), 
                     times = c(length(genes.repem.only),length(overlap))))
dataset2 <- data.frame(value = c(mi_repem_only2, mi_all2), group = factor)
par(mar = c(3, 4.5, 1, 1))
boxplot(value ~ group, dataset2, col = c("antiquewhite1", "lightseagreen"), 
        xlab = NULL, ylab = "Moran's I statistics",outline = FALSE,
        cex.axis = 1.4, cex.lab = 1.6)

##--------------------------------------------------------------------
## Summarized spatial expression patterns
##--------------------------------------------------------------------
library(amap)
source("./funcs/PlotFuncs.R")

counts <- cbind(counts1[overlap,], counts8[overlap,])
location8.trans <- location8
location8.trans$x <- location8$x + max(location1$x)-2
location <- rbind(location1, location8.trans)

vst_count <- var_stabilize(counts) # R function in funcs.R
sig_vst_count <- vst_count[genes_rep_repem, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(apply(counts, 2, sum))))

hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 3
memb <- cutree(hc, k = numC)

# The mean residuals of the three patterns for each location
cent <- NULL
for (k in 1:numC) {
  cent <- cbind(cent, colMeans(sig_vst_res[memb == k, , drop = FALSE]))
}

position_cord <- location
rownames(position_cord) <- rownames(cent)

# The relative residuals
rel_cent <- t(apply(cent, 1, relative_func))
pd <- setNames(cbind.data.frame(position_cord, rel_cent), c("x", "y", paste0("Pattern ", c("III", "II", "I"))))
MBP <- lapply(1:numC, function(x) {
  pattern_plot2(pd, x, xy = T, main = T, titlesize = 2)
})

grid.arrange(grobs = MBP[numC:1], nrow = 3) # 8.5 * 4.75 inches

##--------------------------------------------------------------------
## Spatial expression patterns of randomly selected SVGs
##--------------------------------------------------------------------
source("funcs/PlotFuncs.R")
genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh) 
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap]

randi <- sample(length(genes.repem.only), 20, replace = FALSE)
gene_plot <- genes.repem.only[randi]

# gene_plot <- c("Slc6a13", "Cbln4", "Clip4") # three representative SVGs only identified by STAREG

# Rep1
vst_ct <- var_stabilize(counts1) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)

pltdat <- cbind.data.frame(location1[,1:2],rel_vst_ct)
genetitle <- gene_plot

pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot2(pltdat,x,main=T,titlesize=2.5,title=genetitle[x])})
grid.arrange(grobs=pp, nrow=2)

# Rep8
vst_ct <- var_stabilize(counts8) # R function in funcs.R
sig_vst_ct <- vst_ct[gene_plot, ]
rel_vst_ct <- apply(sig_vst_ct, 1, relative_func)

pltdat <- cbind.data.frame(location8[,1:2],rel_vst_ct)
genetitle <- gene_plot

pp <- lapply(1:(ncol(pltdat)-2),
             function(x){pattern_plot2(pltdat,x,main=T,titlesize=2,title=genetitle[x])})
grid.arrange(grobs=pp, nrow=3)

##--------------------------------------------------------------------
## UMAP clusters
##--------------------------------------------------------------------
library(umap)
library(amap)
source("./R/PlotFuncs.R")

# anscombe variance stabilizing transformation: NB
vst_count <- var_stabilize(counts1) # R function in funcs.R
sig_vst_count <- vst_count[genes.repem.only, ]
sig_vst_res <- t(apply(sig_vst_count, 1, LMReg, T = log(apply(counts1, 2, sum))))

hc <- hcluster(sig_vst_res, method = "euc", link = "ward", nbproc = 1, 
               doubleprecision = TRUE)
numC <- 3
memb <- cutree(hc, k = numC)

sig_counts <- counts1[genes.repem.only,]
scaled_counts <- t(scale(t(sig_counts)))
umap_results <- umap::umap(scaled_counts)

labels <- memb
labels[which(labels == 1)] <- "I"
labels[which(labels == 2)] <- "II"
labels[which(labels == 3)] <- "III"

df <- data.frame(umap_results$layout[,1], umap_results$layout[,2], labels)
# df <- df[which(memb<=4 & memb>=2),]
colnames(df) <- c("UMAP1","UMAP2", "Cluster")
ggplot(df, aes(x = UMAP1, y = UMAP2, color = Cluster), color = labels) + 
  geom_point() + 
  scale_colour_manual(name="Cluster",  
                      values = c("I"="#3cc08f", "II"="#F4A460", "III"="Steelblue"))

##--------------------------------------------------------------------
## Validation gene sets
##--------------------------------------------------------------------
library(rjson)
library(stringr)
bh.only    <- genes_rep_bh[!genes_rep_bh%in%genes_rep_maxp]
repem.only <- genes_rep_repem[!genes_rep_repem%in%genes_rep_maxp]

## The highlighted marker genes in the original study
spatial.genes <- c("Doc2g", "Slc17a7", "Reln", "Cdhr1", "Sv2b", "Shisa3", "Plcxd2", "Nmb", "Uchl1", "Rcan2")
length(intersect(genes_rep_maxp, spatial.genes)) # 7/559
length(intersect(genes_rep_bh, spatial.genes)) # 7/618
length(intersect(genes_rep_repem, spatial.genes)) # 9/1175

length(intersect(bh.only, spatial.genes)) # 0/59
length(intersect(repem.only, spatial.genes)) # 2/616

## Genes related to the main olfactory bulb in Harmonizome database: Allen Brain Atlas
library(rjson)
library(stringr)
Harmonizome_Allen <- fromJSON(paste(readLines("./validation/mouse olfactory bulb/Harmonizome_Allen.json")))
Harmonizome_Allen <- Harmonizome_Allen[["associations"]]
Allen.genes <- NULL
for (i in 1:length(Harmonizome_Allen)){
  Allen.genes   <- c(Allen.genes, Harmonizome_Allen[[i]][["gene"]][["symbol"]])
}
Allen.genes <- Allen.genes[!duplicated(Allen.genes)]
Allen.genes <- tolower(Allen.genes)
Allen.genes <- str_to_title(Allen.genes) # 1894

length(intersect(genes_rep_maxp, Allen.genes)) # 125/559
length(intersect(genes_rep_bh, Allen.genes)) # 132/618
length(intersect(genes_rep_repem, Allen.genes)) # 216/1175

length(intersect(bh.only, Allen.genes)) # 7/59
length(intersect(repem.only, Allen.genes)) # 91/616


## Genes relative to olfactorybulb: BioGPS
Harmonizome_BioGPS <- fromJSON(paste(readLines("./validation/mouse olfactory bulb/Harmonizome_BioGPS.json")))
Harmonizome_BioGPS <- Harmonizome_BioGPS[["associations"]]

BioGPS.genes <- NULL
BioGPS.pvalue <- NULL
for (i in 1:length(Harmonizome_BioGPS)){
  BioGPS.genes <- c(BioGPS.genes, Harmonizome_BioGPS[[i]][["gene"]][["symbol"]])
  BioGPS.pvalue <- c(BioGPS.pvalue, 10^(-abs(Harmonizome_BioGPS[[i]][["standardizedValue"]])))
}
names(BioGPS.pvalue) <- BioGPS.genes
BioGPS.genes <- BioGPS.genes[!duplicated(BioGPS.genes)]
BioGPS.genes <- tolower(BioGPS.genes)
BioGPS.genes <- str_to_title(BioGPS.genes) # 2031

length(intersect(genes_rep_maxp, BioGPS.genes)) # 191/559
length(intersect(genes_rep_bh, BioGPS.genes)) # 203/618
length(intersect(genes_rep_repem, BioGPS.genes)) # 315/1175

length(intersect(bh.only, BioGPS.genes)) # 12/59
length(intersect(repem.only, BioGPS.genes)) # 124/616

##--------------------------------------------------------------------
## GO enrichment analysis
##--------------------------------------------------------------------
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(ggplot2)
library(ggrepel)
# options(connectionObserver = NULL) # run if library(org.Mm.eg.db) failed

genes.all <- bitr(overlap, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh)
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
sum(go.repem.only$p.adjust<0.05) # 245
# sig.go <- filter(go.repem.only, p.adjust<.05)
# dotplot(sig.go,showCategory = 245, size = NULL,font.size = 10, title = "GO enrichment", split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
# barplot(sig.go,showCategory = 245, size = NULL,font.size = 10, title = "GO enrichment", split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")

write.csv(go.repem.only,file = "output/MOB_go_repem_only.csv",quote = FALSE)

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
  dplyr::select(-len) %>%
  
  # Add this info to the initial dataset
  left_join(results, ., by=c("ONTOLOGY"="ONTOLOGY")) %>%
  
  # Add a cumulative position of group
  arrange(ONTOLOGY, BP) %>%
  mutate(BPcum=BP+tot)

axisdf = don %>% group_by(ONTOLOGY) %>% summarize(center=(max(BPcum) + min(BPcum))/2)
annotated = c("extracellular space", "trans-synaptic signaling", "chemical synaptic transmission",
              "collagen-containing extracellular matrix", "cell adhesion molecule binding", "sensory perception",
              "dendrite", "neuronal cell body", "presynapse", "dopamine receptor binding")

ggplot(don, aes(x = BPcum, y=-log10(pvalue))) +
  geom_point(aes(color = as.factor(ONTOLOGY), size = Count), alpha=0.8) +
  scale_colour_manual(name="",  
                      values = c("BP"="Tan", "CC"="DarkSeaGreen", "MF"="#6496D2")) +
  scale_x_continuous(label = axisdf$ONTOLOGY, breaks= axisdf$center) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 13.5)) +     # remove space between plot area and x axis
  geom_hline(yintercept = -log10(0.0033), color = '#545454', size= 1.2, linetype = "dashed") + 
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

genes.overlap    <- intersect(genes_rep_repem, genes_rep_bh) 
genes.repem.only  <- genes_rep_repem[!genes_rep_repem%in%genes.overlap]
genes_repem_only <- bitr(genes.repem.only, fromType = "SYMBOL", toType = c("ENTREZID"), OrgDb = org.Mm.eg.db)
kegg.repem.only <- enrichKEGG(gene         = genes_repem_only$ENTREZID,
                              organism     = 'mmu',
                              universe     = genes.all$ENTREZID,
                              pvalueCutoff = 0.9,
                              qvalueCutoff =0.9)
sum(kegg.repem.only$p.adjust<0.05) # 11

sig.kegg <- filter(kegg.repem.only, p.adjust<.05)
dim(sig.kegg)
kegg.bar <- barplot(sig.kegg,showCategory=20,color = "pvalue")
par(mar = c(4, 4.5, 0.5, 1))
kegg.dot <- dotplot(sig.kegg,showCategory=20,color = "pvalue") + 
  theme( 
    legend.position = c(0.98,0.02),
    legend.justification = c(0.98,0.02),
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 12), # face = "bold",
    axis.title = element_text(size = 15),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank() 
  ) + scale_color_continuous()
