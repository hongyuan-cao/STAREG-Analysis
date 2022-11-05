library(SPARK)
source("./funcs/spark.test.R")

##---------------------------------------------------------------------------------------
## load the 10x human breast cancer data (Section 1)
## https://drive.google.com/drive/folders/1FuvWrHenJuXApEet9Rf0HSuLvL5_LnC8?usp=sharing
##---------------------------------------------------------------------------------------
sp_count1 <- Matrix::readMM("./real_data/10X Visium/Section_1/filtered_feature_bc_matrix/matrix.mtx")
genes <- read.delim("./real_data/10X Visium/Section_1/filtered_feature_bc_matrix/features.tsv", header=FALSE)
barcodes <- read.table("./real_data/10X Visium/Section_1/filtered_feature_bc_matrix/barcodes.tsv", quote="\"", comment.char="")
location1 <- read.csv("./real_data/10X Visium/Section_1/spatial/tissue_positions_list.csv", header=FALSE)

barcodes <- as.vector(t(barcodes))
rownames(sp_count1) <- genes$V2
rownames(location1) <- location1$V1
location1 <- location1[barcodes,5:6]
colnames(location1) <- c("x", "y")
colnames(sp_count1) <- paste(location1$x, location1$y, sep = "x")
rownames(location1) <- colnames(sp_count1)

# SPARK analysis
spark1 <- CreateSPARKObject(counts = sp_count1, location = location1, 
                            percentage = 0.1, min_total_counts = 10)
spark1@lib_size <- apply(spark1@counts, 2, sum)
spark1 <- spark.vc(spark1, covariates = NULL, lib_size = spark1@lib_size, 
                   num_core = 5, verbose = T, fit.maxiter = 500)
spark1 <- spark.test(spark1, check_positive = T, verbose = T)

##-----------------------------------------------------------------------
## load the ST human breast cancer data
##-----------------------------------------------------------------------
counts2 <- read.table("./real_data/Layer2_BC_count_matrix-1.tsv", check.names = F)
location2 <- cbind.data.frame(x = as.numeric(sapply(strsplit(rownames(counts2), split = "x"), "[", 1)), 
                              y = as.numeric(sapply(strsplit(rownames(counts2), split = "x"), "[", 2)))
rownames(location2) <- rownames(counts2)
counts2 = t(counts2)
spark2 <- CreateSPARKObject(counts = counts2, location = location2, 
                            percentage = 0.1, min_total_counts = 10)
spark2@lib_size <- apply(spark2@counts, 2, sum)
spark2 <- spark.vc(spark2, covariates = NULL, lib_size = spark2@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark2 <- spark.test(spark2, check_positive = T, verbose = T)

##-------------------------------------------------------------------------------------
## https://drive.google.com/file/d/1QeMJThlvKpqcFBd4Zm9Quuu4VBomdgjV/view?usp=sharing
##-------------------------------------------------------------------------------------
save(spark1, spark2, file = "./output/HBC.RData")

