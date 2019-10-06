


#' Impute data using the MAGIC algorithm
#' 
#' This is an R wrapper for the orginal python implementation of the Marcov Affinity-based Graph Imputation of Cells (MAGIC) algorithm as described in 
#' Van Dijk, David et al. The default behavior of the orginal python implementation is slightly different and sometimes preferable to the existing R 
#' implementation of MAGIC (Rmagic::magic). magicBatch also includes support for imputation across multiple batches: the optional aff_mat_input argument
#' allows the user to pass any low-dimensional representation of the data, including batch-corrected data, directly to the calculation of the powered marcov 
#' affinity matrix, resulting in cross-batch imputation that preserves batch-specific differences.
#' 
#' @param data An expression matrix where cells correspond to rows and genes correspond to columns
#' @param aff_mat_input A matrix where cells correspond to rows and components correspond to columns. If left unspecified, the affinity matrix calculation is initialized with PCA as in the original implementation.
#' @param pca An integer specifying the number of PCA components that should be used
#' @param t The power to which the marcov matrix is to be raised
#' @param k The number of nearest neighbors used to construct the knn graph
#' @param ka This controls the standard deviation used in the Gaussian kernel width for a given cell, which is set to the distance to the ka-th nearest neighbor.
#' @param epsilon Epsilon parameter used in MAGIC
#' @param rescale Percentile to rescale data to after imputation
#' @param n_diffusion_components Number of diffusion map components to compute.
#' @param python_command A character string passed to the "command" arugment of the system2 function in order to invoke python.
#' @export

magicBatch <- function(data,
                       aff_mat_input = NULL,
                       pca = 20,
                       t = 2,
                       n_diffusion_components = 10,
                       k = 9, 
                       ka = 3, 
                       epsilon = 1, 
                       rescale = 0,
                       python_command = "python3") {
  
  on.exit({print(paste("removing temporary files"))
    file.remove(to.remove)
    print(paste("done") )})
  path2MagScript <- system.file("exec", "MAGIC.py", package="magicBatch")
  path2MagInputData <- tempfile("magicData", fileext = ".csv")
  message("Input Data Path: ", path2MagInputData)
  data.table::fwrite(as.data.frame(data), file = path2MagInputData, row.names = FALSE)
  path2MagOutput <- tempfile("magicOut", fileext = ".csv")
  path2MagOutput2 <- tempfile("magicOut2", fileext = ".csv")
  path2MagOutput3 <- tempfile("magicOut3", fileext = ".csv")

  message("Output Data Paths: ", "\n", path2MagOutput, "\n", path2MagOutput2, "\n", path2MagOutput3)
  to.remove <- c(path2MagInputData, path2MagOutput, path2MagOutput2, path2MagOutput3)
  magParams <- paste("-d", path2MagInputData,
                     "-o", path2MagOutput,
                     "-q", path2MagOutput2,
                     "-s", path2MagOutput3,
                     "-c", n_diffusion_components,
                     "-p", pca,
                     "-t", t,
                     "-k", k,
                     "-ka", ka,
                     "-e", epsilon,
                     "-r", rescale,
                     "-n")
  if(!is.null(aff_mat_input)) {
    path2MagInputAffMat <- tempfile("magicAffMatInput", fileext = ".csv")
    data.table::fwrite(as.data.frame(aff_mat_input), file = path2MagInputAffMat, row.names = TRUE)
    magParams <- paste(magParams, "-a", path2MagInputAffMat)
    to.remove <- c(to.remove, path2MagInputAffMat)
  }
  output <- system2(python_command, args= paste(path2MagScript, magParams, "csv"), stdout=TRUE)
  print(paste(output))
  imputed_data <- as.matrix(data.table::fread(path2MagOutput, header = TRUE)[,-1])
  dimnames(imputed_data) <- dimnames(data)
  diffusion_map <- as.matrix(data.table::fread(path2MagOutput2, header = TRUE)[,-1])
  rownames(diffusion_map) <- rownames(data)
  affinity_matrix <- as.matrix(data.table::fread(path2MagOutput3, header = TRUE)[,-1])
  dimnames(affinity_matrix) <- list(rownames(data), rownames(data))
  return(list(imputed_data = imputed_data,
              diffusion_map = diffusion_map,
              affinity_matrix = affinity_matrix))
}

