########################################################################################
#               Real data guided Simulation function (mouse olfactory bulb)            #
#              Study 1: based on Replicate 1 + Study 2: based on Replicate 8           #
########################################################################################
#' @param tau1: a numeric value, the variance of non-spatial residual error in study 1, 0.2 as default.
#' @param tau2: a numeric value, the variance of non-spatial residual error in study 2, 0.2 as default.
#' @param ipt: a string in {"Pattern I", "Pattern II", "Pattern III"}, indicating the spatial pattern used for simulation. "Pattern II" as default.
#' @param sig1.str: a string, strength of SVG signals for study 1. "Strong" as default, also support "Weak" and "Strong".
#' @param sig2.str: a string, strength of SVG signals for study 2. "Strong" as default, also support "Weak" and "Strong".
#' @param truth1: the true hidden states of study 1.
#' @param truth2: the true hidden states of study 2.
#'
#' @return a list with the elements 
#' \item{pvals1}{the simulated p-values for study 1.} 
#' \item{pvals2}{the simulated p-values for study 2.} 
#' \item{truth1}{the true hidden states of study 1.} 
#' \item{truth2}{the true hidden states of study 2.}


SimuData <-
  function(tau1 = 0.2,
           tau2 = 0.2,
           ipt = "Pattern II",
           sig1.str = "Moderate",
           sig2.str = "Moderate",
           truth1, 
           truth2) {
    require(SPARK)
    source("./funcs/spark.test.R")
    n.gene = length(truth1)
    
    load("./output/MOB_pattern.RData") # obtained from https://drive.google.com/file/d/1l-MBOF3gSCWMi9g99tsDkMGcY0bcrYA4/view?usp=sharing
    
    ## Generate the Study 1 spatial count data based on the SPARK results from the real
    ## mouse olfactory bulb data
    beta1 <- sapply(1:length(spark1@res_vc), function(x) {
      spark1@res_vc[[x]]$coefficients
    })
    
    # -10.46, -9.77, -9.36, -9.08
    nb1 <- sapply(1:4, function(x) {
      log(x * exp(median(beta1)))
    })  
    
    info1 <- cbind(spark1@location, total_counts = spark1@lib_size)
    newN1 <- info1$total_counts
    pattern1 <- datalist1[,ipt]
    grp1 <- as.numeric(pattern1 > mean(pattern1)) + 1
    
    if (sig1.str == "Weak")
      ifc1 <- 2
    if (sig1.str == "Moderate")
      ifc1 <- 3
    if (sig1.str == "Strong")
      ifc1 <- 4
    uu1  <- c(nb1[1], nb1[ifc1], nb1[1])
    
    lambda1 <- sapply(1:n.gene, function(x) {
      truth1[x] * exp(uu1[grp1] + rnorm(length(uu1[grp1]), 0, tau1)) + (!truth1[x]) *
        exp(uu1[3] + rnorm(length(uu1[grp1]), 0, tau1))
    })
    newCt1 <- lapply(1:(n.gene), function(x) {
      rpois(length(lambda1[, x]), newN1 * lambda1[, x])
    })
    countdata1 <- data.frame(do.call(rbind, newCt1))
    rownames(countdata1) <- paste0("gene", 1:nrow(countdata1))
    colnames(countdata1) <- rownames(info1)
    
    ## Compute primary p-values with SPARK
    spark1 <-
      CreateSPARKObject(
        counts = countdata1,
        location = info1[,1:2],
        percentage = 0.1,
        min_total_counts = 10
      )
    
    spark1@lib_size <- info1$total_counts
    spark1 <-
      spark.vc(
        spark1,
        covariates = NULL,
        lib_size = spark1@lib_size,
        num_core = 5,
        verbose = T,
        fit.maxiter = 500
      )
    spark1 <- spark.test(spark1, check_positive = T, verbose = T)
    pvalue1 <- spark1@res_mtest[, "combined_pvalue"]
    closeAllConnections()
    
    ## Generate the Study 2 spatial count data based on the SPARK results from the real
    ## mouse olfactory bulb data
    beta2 <- sapply(1:length(spark8@res_vc), function(x) {
      spark8@res_vc[[x]]$coefficients
    })
    
    # -9.98, -9.29, -8.88, -8.59
    nb2 <- sapply(1:4, function(x) {
      log(x * exp(median(beta2)))
    })

    info2 <- cbind(spark8@location, total_counts = spark8@lib_size)
    newN2 <- info2$total_counts
    pattern2 <- datalist2[,ipt]
    grp2 <- as.numeric(pattern2 > mean(pattern2)) + 1
    
    if (sig2.str == "Weak")
      ifc2 <- 2
    if (sig2.str == "Moderate")
      ifc2 <- 3
    if (sig2.str == "Strong")
      ifc2 <- 4
    uu2  <- c(nb2[1], nb2[ifc2], nb2[1])
    
    lambda2 <- sapply(1:n.gene, function(x) {
      truth2[x] * exp(uu2[grp2] + rnorm(length(uu2[grp2]), 0, tau2)) + (!truth2[x]) *
        exp(uu2[3] + rnorm(length(uu2[grp2]), 0, tau2))
    })
    newCt2 <- lapply(1:(n.gene), function(x) {
      rpois(length(lambda2[, x]), newN2 * lambda2[, x])
    })
    countdata2 <- data.frame(do.call(rbind, newCt2))
    rownames(countdata2) <- paste0("gene", 1:nrow(countdata2))
    colnames(countdata2) <- rownames(info2)
    
    ## Compute primary p-values with SPARK
    spark2 <-
      CreateSPARKObject(
        counts = countdata2,
        location = info2[, 1:2],
        percentage = 0.1,
        min_total_counts = 10
      )
    
    spark2@lib_size <- info2$total_counts
    spark2 <-
      spark.vc(
        spark2,
        covariates = NULL,
        lib_size = spark2@lib_size,
        num_core = 5,
        verbose = T,
        fit.maxiter = 500
      )
    spark2 <- spark.test(spark2, check_positive = T, verbose = T)
    pvalue2 <- spark2@res_mtest[, "combined_pvalue"]
    closeAllConnections()
    
    return(list(
      pvals1 = pvalue1,
      pvals2 = pvalue2,
      truth1 = truth1,
      truth2 = truth2))
  }