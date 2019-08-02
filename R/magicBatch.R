



dat.comb <- readRDS("~/git/magic/data/dat.comb.rds")


magicOg <- function(data,
                    aff_mat_input = reducedDim(dat.comb, "MNN.cc"),
                    mag_pca = 20,
                    mag_t = 2,
                    mag_k = 9, 
                    mag_ka = 3, 
                    mag_epsilon = 1, 
                    mag_rescale = 99,
                    python_command = "python3",
                    path2MagScript = "~/git/magic/src/magic/MAGIC.py") {
  
  path2MagInputData <- tempfile("magicData", fileext = ".csv")
  data.table::fwrite(as.data.frame(data), file = path2MagInputData, row.names = T)
  path2MagOutput <- tempfile("magicOut", fileext = ".csv")
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
    data.table::fwrite(as.data.frame(aff_mat_input), file = path2MagInputAffMat, row.names = T)
    magParams <- paste(magParams, "-a", path2MagInputAffMat)
    to.remove <- c(to.remove, path2MagInputAffMat)
  }
  output <- system2(python_command, args= c(path2MagScript, paste(magParams, "csv")), stdout=TRUE)
  print(paste(output))
  imputedData <- as.matrix(data.table::fread(path2MagOutput, header = T)[,-1])
  imputedData[is.na(imputedData)] <- 0
  dimnames(imputedData) <- dimnames(data)
  on.exit(file.remove(to.remove))
  return(imputedData)
}

test <- magicOg(data = assay(dat.comb, "logcounts"),
        aff_mat_input = reducedDim(dat.comb, "MNN.cc"))
