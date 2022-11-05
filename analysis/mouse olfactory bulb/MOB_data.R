###############################################################
#             Mouse olfactory bulb data analysis              #
###############################################################
library(SPARK)

##------------------------------------------------------------------------
## Load the spatial MOB data Replicate 1 and Replicate 8 and apply SPARK 
##------------------------------------------------------------------------
# read the raw counts (spatial data)
counts1 <- read.table("./real_data/Rep1_MOB_count_matrix-1.tsv", check.names = F)
counts8 <- read.table("./real_data/Rep8_MOB_count_matrix-1.tsv", check.names = F)
counts1 <- t(counts1)
counts8 <- t(counts8)

# Replicate 1
location1 <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(counts1), split = "x"), "[", 1)), 
                              y = as.numeric(sapply(strsplit(colnames(counts1), split = "x"), "[", 2)))
rownames(location1) <- colnames(counts1)

spark1 <- CreateSPARKObject(counts = counts1, location = location1, 
                            percentage = 0.1, min_total_counts = 10)
spark1@lib_size <- apply(spark1@counts, 2, sum)
spark1 <- spark.vc(spark1, covariates = NULL, lib_size = spark1@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark1 <- spark.test(spark1, check_positive = T, verbose = T)

# Replicate 8
location8 <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(counts8), split = "x"), "[", 1)), 
                              y = as.numeric(sapply(strsplit(colnames(counts8), split = "x"), "[", 2)))
rownames(location8) <- colnames(counts8)
spark8 <- CreateSPARKObject(counts = counts8, location = location8, 
                            percentage = 0.1, min_total_counts = 10)
spark8@lib_size <- apply(spark8@counts, 2, sum)
spark8 <- spark.vc(spark8, covariates = NULL, lib_size = spark8@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark8 <- spark.test(spark8, check_positive = T, verbose = T)


##----------------------------------------------------------------------------
## https://drive.google.com/file/d/1kUr8wW1NheyevfxzP7Fx6Wy3Nndu6j51/view?usp=sharing
save(counts1, counts8, location1, location8, spark1, spark8, file = "./output/MOB.RData")
