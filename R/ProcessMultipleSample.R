#' @@ ProcessMultipleSample
#' @description Several steps to deal with users uploaded multi-sample expression profile:
#' calculating and normalizing ES score of samples, combine the normalized ES score from
#' positive set and negative set of genes belong to the same step by calculating the difference between them,
#' then change the format of results for visualization on web.
#'
#' @param expression.from.users A data.frame contains expression data from users.
#' @param signatureList List of signature genes, default by 'signature.GeneSymbol.list'.
#' @param perm.times Numeric value indicating times of permutation, default by 100.
#' @param signature.annotation A data.frame contains step and function direction information about signature genes,
#' default by 'annotation'.
#' @param type.of.data Character indicating source of expression data, 'Microarray'or 'RNA-seq'.
#' @param format.of.file Character indicating format of RNA-seq expression data, 'TPM' or 'Count'.
#' @param gene.length.path Character represting the storage path of the gene length file.

ProcessMultipleSample <- function(expression.from.users,
                                  signatureList = signature.GeneSymbol.list,
                                  perm.times = 100,
                                  signature.annotation = annotation,
                                  type.of.data = c("Microarray", "RNA-seq"),
                                  format.of.file = c("TPM","Count"),
                                  gene.length.path){
  ##########
  ##---------1. signature match(at least 1----------------------------------------------------------------------
  ##########
  #1.intersect number of every signature sets
  separate.remain <- sapply(signatureList, function(i){
    num <- length(intersect(i, rownames(expression.from.users)))
    return(num)
  })

  show.num <- paste0("(", separate.remain, "/", sapply(signatureList, length), ")")
  names(show.num) <- names(signatureList)

  #2.combine the positive genes and negative genes
  positive.l <- grep("positive", names(signatureList))
  #mix.show.num is for visualization
  mix.show.num <- data.frame(separate.remain = separate.remain, signature.num = sapply(signatureList, length),
                             step = c(1,2,2,3,3,4:20,21,21,22,22,23,23))


  l0 <- which(separate.remain == 0)
  #3.if one direction of a step matched 0 genes, the other direction will be also discarded
  a <- which(l0 %in% positive.l)
  b <- which(l0 %in% (positive.l-1))
  if(length(a > 0) | length(b > 0)){
    delete.signature <- c(l0[-c(a,b)],l0[a],l0[a]-1,l0[b],l0[b]+1)
  }else{
    delete.signature <- l0
  }

  if(length(delete.signature) > 0){
    signatureList <- signatureList[-delete.signature]

    mix.show.num <- mix.show.num[-delete.signature,]

    get.separate <- unlist(strsplit(names(delete.signature),"\\."))
    delete.step <- intersect(get.separate[grep("Step", get.separate)], c("Step2","Step3","Step5","Step6","Step7"))
    delete.step <- gsub("Step","",delete.step)
    #4.delete needless siganture sets
    if(length(delete.step) > 0){
      signature.annotation <- signature.annotation[-which(signature.annotation[,4] %in% delete.step), ]
    }

  }
  print(delete.signature)
  print(show.num)

  mix.show.num$step <- as.factor(mix.show.num$step)
  mix.show.num <- paste0("(", tapply(mix.show.num[,1], mix.show.num$step, sum), "/",
                         tapply(mix.show.num[,2], mix.show.num$step, sum), ")")
  positive.l <- grep("positive", names(signatureList))
  if(length(positive.l) > 0){
    names(mix.show.num) <- names(signatureList)[-positive.l]
    names(mix.show.num) <- gsub(".negative", "", names(mix.show.num))
  }else{
    names(mix.show.num) <- names(signatureList)
  }

  mix.show.num <- paste(names(mix.show.num), mix.show.num, sep = " ")



  ##########
  ##---------2. ssGSEA process---------------------------------------------------------------------------------
  ##########
  #1. Data preprocessing
  #log2 transformation of expression data, and plus 1 for all values to make accurate rank later
  if(type.of.data == "RNA-seq" & format.of.file == "TPM"){
    normalized.expression.from.users <- as.matrix(log2(expression.from.users+1)+1)

  }else if(type.of.data == "RNA-seq" &  format.of.file == "Count"){
    #transform count to tpm
    normalized.expression.from.users <- count2Tpm(countName = expression.from.users, filePath = gene.length.path, sample="multiple")
    print(paste("count to TPM down---dim(normalized.expression.from.users)", dim(normalized.expression.from.users)))
    normalized.expression.from.users <- as.matrix(log2(normalized.expression.from.users+1)+1)

  }else{  #if the microarray data was log-transformed
    if(max(expression.from.users) < 50){
      normalized.expression.from.users <- as.matrix(expression.from.users+1)
    }else{
      normalized.expression.from.users <- as.matrix(log2(expression.from.users+1)+1)
    }

  }

  #2. build the random matrix 'permutation.exp', permutation of matrix 'expression.from.users'
  sample.number <- ncol(expression.from.users)
  permutation.exp <- matrix(rep(normalized.expression.from.users, times = (perm.times+1)),
                            nrow = nrow(normalized.expression.from.users))
  #3.ssGSEA for the whole matrix 'permutation.exp'
  for(i in 1:perm.times){
    a <- ssGSEA.MultipleSamples.Permutation(exp.profile = normalized.expression.from.users)
    permutation.exp[, i*sample.number + 1:sample.number] <- a
  }
  rownames(permutation.exp) <- rownames(normalized.expression.from.users)

  exp.n <- ncol(permutation.exp)
  exp.m <- nrow(permutation.exp)
  #(1)get the superimposed matrix
  superposition.rank.matrix <- superpositionRank(raw.matrix = permutation.exp, m = exp.m, n = exp.n)
  #(2)get rank by column of matrix 'permutation.exp'
  permutation.exp.rank <- rankMatrixByCol(superposition.rank = superposition.rank.matrix, m = exp.m, n = exp.n)

  #activity score of signature sets one by one
  permutation.score <- t(sapply(signatureList, function(sig){
    gSetIdx <- which(rownames(permutation.exp) %in% sig)
    #(3)vector of '0'and '1' show the position of signature genes in 'permutation.exp'
    t1 <- numeric(exp.m)
    t1[gSetIdx] <- 1
    geneSet.position <- matrix(rep(t1, times = exp.n), ncol = exp.n)
    #(4)activity score of one signature set
    score.for.oneset <- oneSetssGSEA.ES(exp.rank = permutation.exp.rank, superposition.rank.es = superposition.rank.matrix,
                                        gSet.position.matrix = geneSet.position, alpha=0.25, m.num = exp.m, n.num = exp.n)

    return(score.for.oneset)
  }))
  print("permutation.score over!")



  ##########
  ##---------3.normalization of activity score-----------------------------------------------------------------------------
  ##########
  #1. whether activity score of real sample and permutation samples(N+1，2N+1 ... perm.times*N+1) are of the same sign
  #if different, set the score of permutation samples to 'NA'
  for(t in 1:perm.times){
    permutation.oneExp.score <- permutation.score[, t*sample.number + 1:sample.number]
    r <- which(permutation.score[,1:sample.number] * permutation.oneExp.score < 0)
    permutation.oneExp.score[r] <- NA
    permutation.score[, t*sample.number + 1:sample.number] <- permutation.oneExp.score
  }

  #2. z-score normalization
  zScore.matrix <- sapply(1:sample.number, function(s){
    #activity score of real sample and permutation samples(N+1，2N+1 ... perm.times*N+1)
    perm.matrix.forOneSample <- permutation.score[,((0:perm.times) * sample.number + s)]
    #z-score normalization
    zScore.forOneSample <- apply(perm.matrix.forOneSample, 1, function(x) {
      if(length(which(!is.na(x))) < 2){
        zScore <- x[1]
      }else{

        zScore <- scale(x)[1]
      }
      return(zScore)
    })
    return(zScore.forOneSample)
  })
  print("zScore over!")

  #3. step 2,3,5,6,7 z-score positive - negative
  positive.l <- grep("positive", rownames(zScore.matrix))
  if(length(positive.l) > 0){
    #change the sign of score of negative genes
    zScore.matrix[positive.l-1, ] <- zScore.matrix[positive.l-1, ]*(-1)
    direction <- 1:length(signatureList)
    direction[positive.l] <- direction[positive.l-1]
  }else{
    direction <- 1:length(signatureList)
  }
  print("direction over!")

  #rowsum to conbine positive genes score and negative genes score
  ssGSEA.normalized.score <- rowsum(zScore.matrix, direction)
  if(length(positive.l) > 0){
    rownames(ssGSEA.normalized.score) <- rownames(zScore.matrix)[-positive.l]
    rownames(ssGSEA.normalized.score) <- gsub(".negative", "", rownames(ssGSEA.normalized.score))
  }else{
    rownames(ssGSEA.normalized.score) <- rownames(zScore.matrix)
  }
  colnames(ssGSEA.normalized.score) <- colnames(expression.from.users)

  ssGSEA.normalized.score <- round(ssGSEA.normalized.score, 3)

  zScore.ssGSEA.normalized.score <- t(apply(ssGSEA.normalized.score, 1, scale))
  colnames(zScore.ssGSEA.normalized.score) <- colnames(expression.from.users)

  #expression of signature expression---------------------------------------------------------------------------------------------
  #(1)heatmap
  anno <- unique(signature.annotation[,3:4])
  inter <- intersect(signature.annotation[,3], rownames(expression.from.users))
  print(paste0("intersect number: ", length(inter)))
  anno <- anno[which(anno[,1] %in% inter),]
  print(paste("dim(anno)---", dim(anno)))

  annotation_row_length <- tapply(anno[,2], anno[,2], length)
  print(paste("sum of annotation_row_length: ", sum(annotation_row_length)))

  expression.profile <- as.matrix(expression.from.users[anno[,1], ])
  print(paste("dim(expression.profile): ",dim(expression.profile)))
  storage.mode(expression.profile) <- "numeric"

  heatmap.expression <- expression.profile
  rownames(heatmap.expression) = 1:nrow(heatmap.expression)

  annotation_row = data.frame(Steps = factor(rep(paste0("Step", names(annotation_row_length)), times = annotation_row_length)))
  rownames(annotation_row) = rownames(heatmap.expression)
  ann_colors <- list(Steps = c(Step1 = "#66C2A5", Step2="#FC8D62", Step3="#8DA0CB", Step4="#E78AC3",
                               Step5 = "#A6D854", Step6="#FFD92F", Step7="#E5C494"))
  #filter outliers
  heatmap.expression <- log2(heatmap.expression+1)
  print(summary(as.vector(heatmap.expression)))
  print(which(apply(heatmap.expression, 1, sd) == 0))
  scale.expression <- t(apply(heatmap.expression,1,scale))
  colnames(scale.expression) <- colnames(heatmap.expression)
  scale.expression <- filterOutliers(scale.expression)

  p <- pheatmap::pheatmap(scale.expression, cluster_row = FALSE, cluster_col = TRUE,
                          col=colorRampPalette(c("navy", "white", "firebrick3"))(100),
                          annotation_row = annotation_row,
                          main = "Signature Genes Expression", labels_row = "", labels_col = NULL, scale = "none",
                          annotation_colors = ann_colors, border_color = NA)

  ##########
  ##---------5.return----------------------------------------------------------------------------------------------------

  return(list(separate.remain = separate.remain,
              show.num = show.num,
              delete.signature = delete.signature,
              ssGSEA.normalized.score = ssGSEA.normalized.score,
              zScore.matrix = zScore.matrix,
              zScore.ssGSEA.normalized.score = zScore.ssGSEA.normalized.score,
              pheatmap = list(scale.expression = scale.expression,
                              fig = p,
                              annotation_row = annotation_row,
                              ann_colors <- list(Steps = c(Step1 = "#66C2A5", Step2="#FC8D62", Step3="#8DA0CB", Step4="#E78AC3",
                                                           Step5 = "#A6D854", Step6="#FFD92F", Step7="#E5C494"))),
              signature.expression = expression.profile))
}
