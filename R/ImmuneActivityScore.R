


#' @@ ImmuneActivityScore
#'
#' @description The main function to calculating immune activity levels based on ssGSEA algorithm,
#' both for multiple samples and single sample profile.
#' @param expr TPM matrix derived from scRNAseq,
#' or normalized gene expression matrix from bulk samples;
#' rows are genes, cols are cells or samples.
#' @param sampleNumber Numeric value indicating number of samples in profile.
#' @param permTimes Numeric value indicating times of permutation, default by 100.
#' @param type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'.
#' @param format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'.
#' @note R package 'pheatmap' is required.
#' @returnType
#' @return NULL
#'
#' @author Liwen Xu

ImmuneActivityScore <- function(expr,sampleNumber, permTimes, type.of.data, format.of.file){

  print("set.seed");
  set.seed(1:100);

  if(sampleNumber > 1){
    result <- ProcessMultipleSample(expression.from.users = expr,
                                    signatureList = signature.GeneSymbol.list,
                                    perm.times = permTimes,
                                    signature.annotation = annotation,
                                    type.of.data=type.of.data,
                                    format.of.file=format.of.file)

  }else{
    result <- ProcessSingleSample(expression.from.users = expr,
                                  signatureList = signature.GeneSymbol.list,
                                  perm.times = permTimes,
                                  signature.annotation = annotation,
                                  type.of.data=type.of.data,
                                  format.of.file=format.of.file)
  }
  return(result)

}
