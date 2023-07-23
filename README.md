# STAREG-Analysis
 All scripts performed in our experiments for replicability analysis of spatially variable gene detection.

## Installation

```R
## Install dependency packages if necessary
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("qvalue")
install.packages(c("Iso"))

## Install STAREG
install.packages("devtools")
devtools::install_github("hongyuan-cao/STAREG")

## Load STAREG
library(STAREG)
```

## An analysis example

We illustrate the replicability analysis of SVG detection using two replicates of the SRT data from mouse olfactory bulb and the ``SPARK`` package.

First, load the mouse olfactory bulb data in two replicates, and separately analyze the two datasets with the ``SPARK`` method.

```R
## Load the count data for mouse olfactory bulb in Replicate 1 and Replicate 8 downloaded from the Spatial Research webset at https://www.spatialresearch.org/resources-published-datasets/doi-10-1126science-aaf2403
counts1 <- read.table("./Rep1_MOB_count_matrix-1.tsv", check.names = F)
counts8 <- read.table("./Rep8_MOB_count_matrix-1.tsv", check.names = F)
counts1 <- t(counts1)
counts8 <- t(counts8)

## Analyze Replicate 1 dataset with SPARK (referring to the SPARK pacakge)
library(SPARK)
location1 <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(counts1), split = "x"), "[", 1)), 
                              y = as.numeric(sapply(strsplit(colnames(counts1), split = "x"), "[", 2)))
rownames(location1) <- colnames(counts1)
spark1 <- CreateSPARKObject(counts = counts1, location = location1, 
                            percentage = 0.1, min_total_counts = 10)
spark1@lib_size <- apply(spark1@counts, 2, sum)
spark1 <- spark.vc(spark1, covariates = NULL, lib_size = spark1@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark1 <- spark.test(spark1, check_positive = T, verbose = T)

## Analyze Replicate 8 dataset with SPARK
location8 <- cbind.data.frame(x = as.numeric(sapply(strsplit(colnames(counts8), split = "x"), "[", 1)), 
                              y = as.numeric(sapply(strsplit(colnames(counts8), split = "x"), "[", 2)))
rownames(location8) <- colnames(counts8)
spark8 <- CreateSPARKObject(counts = counts8, location = location8, 
                            percentage = 0.1, min_total_counts = 10)
spark8@lib_size <- apply(spark8@counts, 2, sum)
spark8 <- spark.vc(spark8, covariates = NULL, lib_size = spark8@lib_size, 
                   num_core = 10, verbose = T, fit.maxiter = 500)
spark8 <- spark.test(spark8, check_positive = T, verbose = T)
```

Extract the well-calibrated p-values from the analysis results of the two studies and match them by gene.

```R
p1 <- spark1@res_mtest$combined_pvalue
names(p1) <- rownames(spark1@counts)
p2 <- spark8@res_mtest$combined_pvalue
names(p2) <- rownames(spark8@counts)

overlap <- intersect(names(p1),names(p2))
pvals1 = p1[overlap]
pvals2 = p2[overlap]
```

Make replicability analysis of the two studies using STAREG, and output the replicable SVGs across the two SRT studies from mouse olfactory bulb.

```R
library(STAREG)
alpha <- 0.05
rep.obj <- stareg(pvals1, pvals2)
rep.svgs <- overlap[which(rep.obj$fdr <= alpha)]
```

## Data and reproducibility

All the R functions for the realistic simulations and real data analysis are contained in “funcs” folder.

All the R scripts to reproduce the realistic simulations are deposited in “real_simulation” folder.

The SRT data used in the real data analysis can be downloaded from the links provided in the *Data availability* section in the manuscript, and all the raw count data is provided in the “real_data” folder and [here](https://drive.google.com/drive/folders/1dtQiIJNn4hgay3OIcISvZ_IQb3Cun0ig?usp=sharing). 

Some intermediate results for the real data analysis can be found in “output” folder and [here](https://drive.google.com/drive/folders/11-v02JBVBL4N3KqhVOBtYBBX5iW09VCj?usp=sharing). 

All the R scripts that can be used to reproduce the results for real data analysis are summarized in “analysis” folder. 

Some reference gene sets used for validations of the real data analysis results are deposited in “validation” folder.



