library(SPARK) 

##--------------------------------------------------------------------------------------
## Slide-seq mouse cerebellum data
## https://drive.google.com/drive/folders/1b9pzmzJWUmr6yGpGWON44yC7shLc7cv2?usp=sharing
##--------------------------------------------------------------------------------------
counts1 <- read.csv("./real_data/Slideseq-Puck_180819_12/Slide-seq_counts.csv", header = FALSE)
counts1 <- as.matrix(counts1)
genes1 <- as.vector(read.csv("./real_data/Slideseq-Puck_180819_12/GeneNames.csv", header = FALSE))
location1 <- read.csv("./real_data/Slideseq-Puck_180819_12/BeadLocationsForR.csv", header = TRUE)
clusters1 <- read.csv("./real_data/Slideseq-Puck_180819_12/AnalogizerClusterAssignments.csv", header = TRUE)
rownames(counts1) <- t(genes1) # 19782 * 32701
colnames(location1) <- c("barcodes", "x", "y")
colnames(counts1) <- paste(location1$x, location1$y, sep = "x")
rownames(location1) <- colnames(counts1)
rownames(clusters1) <- rownames(location1[which(location1$barcodes %in% clusters1$Var1),])
counts1 <- counts1[,rownames(clusters1)] # 19782 * 28352
location1 <- cbind(location1[rownames(clusters1),], cluster = clusters1$x)
mt_idx  <- grep("mt-",rownames(counts1)) # 26 mitochondrial genes
counts1 <- counts1[-mt_idx,] 
rm(clusters1, genes1)

# SPARK-X analysis
sparkx1 <- sparkx(counts1, location1[,2:3], numCores=1, option="mixture") # 18082 * 28352

##--------------------------------------------------------------------------------------
## Slide-seqV2 mouse cerebellum data
## https://drive.google.com/file/d/11bqirFOjBubWqM1HFqetj-bgRiLFcTXy/view?usp=sharing
##--------------------------------------------------------------------------------------
load("./real_data/SlideseqV2_ROI.rds")
# remove mt genes
mt_idx   <- grep("mt-",rownames(sp_count))
counts2 <- as.matrix(sp_count[-mt_idx,]) 
location2 = as.data.frame(location)
rownames(location2) <- colnames(counts2)
rm(mt_idx, sp_count, location)

# SPARK-X analysis
sparkx2 <- sparkx(counts2, location2, numCores=1, option="mixture") # 20117 * 11626

##--------------------------------------------------------------------------------------
## https://drive.google.com/file/d/13a-GPtPDiwIBNjVbOOOB64x6fQfAh2hn/view?usp=sharing
save(counts1, counts2, location1, location2, sparkx1, sparkx2, file = "./output/MC.RData")



