## Function spark.test() from SPARK package, with a slight modification to close the connections after each iteration.


spark.test <- function(object, 
                       kernel_mat = NULL, 
                       check_positive = TRUE, 
                       verbose = TRUE) {
  
  ## the variables input
  if(!is.null(kernel_mat) & !is.list(kernel_mat)){ stop("SPARK.TEST::kernel_mat must be a list, please provide list type!")}
  ## testing one kernel at a time
  res_kernel <- list()
  res_pval <- NULL
  if(is.null(kernel_mat)){ #calculate kernel matrix using location information
    # Euclid distance, and compute the range of parameter l
    ED <- as.matrix(dist(object@location[ ,1:2]))
    lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
    ## total 10 kernels, i.e., 5 gaussian and 5 periodic kernels
    for(ikernel in c(1:5) ){
      # Gaussian kernel
      cat(paste0("## testing Gaussian kernel: ",ikernel,"...\n"))
      kernel_mat <- exp(-ED^2/(2*lrang[ikernel]^2))
      object <- spark.test_each(object, kernel_mat=kernel_mat, check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
      rm(kernel_mat)
      
      closeAllConnections()
      # Periodic kernel
      cat(paste0("## testing Periodic kernel: ",ikernel,"...\n"))
      kernel_mat <- cos(2*pi*ED/lrang[ikernel])
      object <- spark.test_each(object, kernel_mat=kernel_mat, check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
      rm(kernel_mat)
      closeAllConnections()
    }# end for ikernel
    colnames(res_pval) <- paste0(c("GSP","COS"), rep(1:5,each=2))
    rownames(res_pval) <- rownames(object@counts)
  }else{# kernel_mat is a list provided by user
    for(ikernel in 1:length(kernel_mat) ){
      # pre-defined kernels
      cat(paste0("## testing pre-defined kernel: ",ikernel,"...\n"))
      object <- spark.test_each(object, kernel_mat=kernel_mat[[ikernel]], check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
    }# end for ikernel
    colnames(res_pval) <- paste0("kernel", 1:length(kernel_mat) )
    rownames(res_pval) <- rownames(object@counts)
  }# end fi
  
  object@res_stest <- res_kernel
  ## integrate ten p-values into one
  num_pval <- ncol(res_pval)
  num_gene <- nrow(res_pval)
  ## compute the weight to p-values integrate
  if(is.null(object@weights)){
    weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
  }else if(!is.matrix(object@weights)){
    weights <- as.matrix(object@weights)
  }else {
    weights <- object@weights
  }# end fi
  combined_pvalue <- CombinePValues(res_pval, weights=weights)
  object@res_mtest <- data.frame(res_pval, combined_pvalue = combined_pvalue,  adjusted_pvalue = p.adjust(combined_pvalue, method="BY") )
  # return the results
  return(object)
}# end function score test
