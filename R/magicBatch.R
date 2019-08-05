
#' Impute data using the MAGIC algorithm
#' 
#' This is an R interface for the orginal python implementation of the Marcov Affinity-based Graph Imputation of Cells (MAGIC) algorithm as described in 
#' Van Dijk, David et al. The default behavior of the orginal python implementation is slightly different and sometimes preferable to the existing R 
#' implementation of MAGIC (Rmagic::magic). magicBatch also includes support for imputation across multiple batches: the optional aff_mat_input argument
#' allows the user to pass any low-dimensional representation of the data, including batch-corrected data, directly to the calculation of the powered marcov 
#' affinity matrix, resulting in cross-batch imputation that preserves batch-specific differences.
#' 
#' @param data An expression matrix where cells correspond to rows and genes correspond to columns
#' @param aff_mat_input A matrix where cells correspond to rows and components correspond to columns. If left unspecified, the affinity matrix calculation is initialized with PCA as in the original implementation.
#' @param mag_pca An integer specifying the number of PCA components that should be used
#' @param mag_t The power to which the marcov matrix is to be raised
#' @param mag_k The number of nearest neighbors used to construct the knn graph
#' @param mag_ka This controls the standard deviation used in the Gaussian kernel width for a given cell, which is set to the distance to the ka-th nearest neighbor.
#' @param mag_epsilon Epsilon parameter used in MAGIC
#' @param mag_rescale Percentile to rescale data to after imputation
#' @param python_command A character string passed to the "command" arugment of the system2 function in order to invoke python.
#' @export

magicBatch <- function(data,
                    aff_mat_input = NULL,
                    mag_pca = 20,
                    mag_t = 2,
                    mag_k = 9, 
                    mag_ka = 3, 
                    mag_epsilon = 1, 
                    mag_rescale = 99,
                    python_command = "python3") {
  
  on.exit({print(paste("removing temporary files"))
          file.remove(to.remove)
          print(paste("done") )})
  path2MagScript <- system.file("exec", "MAGIC.py", package="magicBatch")
  path2MagInputData <- tempfile("magicData", fileext = ".csv")
  print(paste("Input Data Path:", path2MagInputData))
  data.table::fwrite(as.data.frame(data), file = path2MagInputData, row.names = FALSE)
  path2MagOutput <- tempfile("magicOut", fileext = ".csv")
  print(paste("Output Data Path: ", path2MagOutput))
  to.remove <- c(path2MagInputData, path2MagOutput)
  magParams <- paste("-d", path2MagInputData,
                     "-o", path2MagOutput,
                     "-p", mag_pca,
                     "-t", mag_t,
                     "-k", mag_k,
                     "-ka", mag_ka,
                     "-e", mag_epsilon,
                     "-r", mag_rescale,
                     "-n")
  if(!is.null(aff_mat_input)) {
    path2MagInputAffMat <- tempfile("magicAffMatInput", fileext = ".csv")
    data.table::fwrite(as.data.frame(aff_mat_input), file = path2MagInputAffMat, row.names = TRUE)
    magParams <- paste(magParams, "-a", path2MagInputAffMat)
    to.remove <- c(to.remove, path2MagInputAffMat)
  }
  output <- system2(python_command, args= c(path2MagScript, paste(magParams, "csv")), stdout=TRUE)
  print(paste(output))
  imputedData <- as.matrix(data.table::fread(path2MagOutput, header = TRUE)[,-1])
  imputedData[is.na(imputedData)] <- 0
  dimnames(imputedData) <- dimnames(data)
  return(imputedData)
}


