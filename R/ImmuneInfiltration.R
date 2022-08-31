#' @title
#' Immune cell infiltration proportion calculate module.
#'
#' @aliases CIBERSORT
#'
#' @description
#' A reference function to produce the calculate cell infiltration proportions results
#' and interactive picture data of the TIP website.
#'
#' @author DengChunyu, Haerbin Medical University, dengcyelena@gmail.com
#' Requirements:
#'       R v3.0 or later. (dependencies below might not work properly with earlier versions)
#'       install.packages('e1071')
#'       install.packages('parallel')
#'       install.packages('preprocessCore')
#'       If preprocessCore is not available in the repositories you have selected, run the following:
#'           source("http://bioconductor.org/biocLite.R")
#'           biocLite("preprocessCore")
#' Windows users using the R GUI may need to Run as Administrator to install or update packages.
#'
#' @example
#'   Example In R:
#'       source('CIBERSORT.server.R')
#'       library('e1071')
#'       library('parallel')
#'       library('preprocessCore')
#'
#'       CIBERSORT.server(codePath="",
#'       filePath="",
#'       signaturePath="",
#'       saveDir="",
#'       perm = 100,
#'       CHIPorRNASEQ="Microarray",
#'       sample="multiple");
#'
#' @param data data
#' @param perm：
#' The NO. permutation. Default by 100.
#'
#' @param CHIPorRNASEQ:
#' The source of data,'RNAseq' or 'Microarray' to choose. Default by "RNA-seq".
#'
#' @param sample：
#' The number of upload sample, "multiple" or "single" sample for user to choose. Default by "multiple".
#'
#' @param dataType:
#' The data type of upload data, "count" or "TPM" for user to choose. Default by "TPM".
#'

CIBERSORT.server <- function(data, perm = 100,CHIPorRNASEQ="RNA-seq", sample="multiple",dataType="TPM"){

  #choose one type of signature matrix we next to use,RNA-seq data for LM14 matrix;Microarray data for LM22 matrix.

  if(CHIPorRNASEQ=="RNA-seq"){

    sig_matrix = LM14

  }else if(CHIPorRNASEQ=="Microarray"){

    sig_matrix = LM22

  }
   print("The signature matrix is import!");

  # run the function to produce result files.
  Result <- CIBERSORT_server(mixture_file= data,
                             sig_matrix = sig_matrix ,
                             perm=perm,
                             CHIPorRNASEQ=CHIPorRNASEQ,
                             sample=sample)

  print("CIBERSORT.server is end");
}



#' @title
#' The web server file about immune infiltration
#'
#' @description
#' This function is to produce all the results file used to produce interactive graphics by web server.
#' We reference the result of function 'CIBERSORT' to produce graphics files of web server need. Including
#' radar plot,stackplot and pie plot.
#'
#' @author DengChunyu
#'
#' @details
#' The interactive graphics files are used to web to produce beautiful and powerful interactive graphics. The format
#' of these file is specific for graphics.
#'
#' @param mixture_file
#' Input data of user upload.
#'
#' @param sig_matrix
#' Input data of signature matrix.
#'
#' @param perm
#' The NO. permutations,100 by defaut.
#'
#' @param QN
#' TRUE by defaut ,Perform quantile normalization or not.
#'
#' @param CHIPorRNASEQ
#' The mixture_file are "Microarray" data or "RNA-seq" data, for user to choose.
#'
#'
#' @param sample
#' User choose 'single' and 'multiple' samples.
#'


CIBERSORT_server <- function(mixture_file,
                             sig_matrix,
                             perm=100,
                             QN=TRUE,
                             CHIPorRNASEQ="RNA-seq",
                             sample="multiple"){

  #Two processing path for the choice of parameter 'CHIPorRNASEQ'

  if(CHIPorRNASEQ=="Microarray"){

    # The result of infiltration proportion.

    Results<-CIBERSORT(sig_matrix=sig_matrix,
                       mixture_file=mixture_file,
                       perm=perm,QN=T,
                       sigMat="LM22",sample=sample)

    # Change the sort of cell types, for the Aesthetic quality of radar plot.

    cell_type<-colnames(Results)[1:(ncol(Results)-3)]

    sort_name_lm22<-c(16,4,21,5,22,10,2,7,11,20,18,17,13,9,1,8,12,15,3,19,14,6)

    cell_type <- cell_type[sort_name_lm22]

    n1<-paste0('["',rownames(Results)[1],'"')

    # Color distribution of 22 immune cell type
    color_lm22<-c("#D03D00","#C3627F","#C84A51","#ACCC69","#969F39","#B3C44C","#7EA33F","#477034","#77B35C","#57B473","#6A84BA","#6067A8","#4BB7A1","#8F88B8","#D2A6CB","#BAB8CC","#5B9EB7","#A7D1E4","#E9D9C3","#C1DDD9","#4FAC6F","#405C70")

    #n m :radar plot; n1 m1 :stackplot

    #prepare for radar plot
    Results_radar<-Results[,sort_name_lm22]

    n<-"name"

    for(j in 1:22){
      n<- paste(n,cell_type[j],sep = ",")
    }

    #Two processing path for the choice of parameter 'sample'
    if(sample=="multiple"){

      #radar
      for(i in 1:dim(Results_radar)[1]){

        m<-"Sample"
        for(j in 1:22)
        {
          m<- paste(m,Results_radar[i,j],sep = ",")

        }
        #radar output
        z<-matrix(c(n,m),ncol=1)

      }

      #stackplot
      for(i in 2:dim(Results)[1]){
        n1<-paste0(n1,',"',rownames(Results)[i],'"')
      }
      n1<-paste0(n1,"]")
      m1_all<-list()
      for(i in 1:22){
        m1<-paste0('{"name":"',colnames(Results)[i],'","color":"',color_lm22[i],'","data":[',Results[1,i])
        for(j in 2:dim(Results)[1]){
          m1<-paste(m1,Results[j,i],sep = ",")
        }
        m1<-paste0(m1,"]}")
        m1_all[i]<-m1
      }
      #stackplot output
      m1_all<-paste0(unlist(m1_all),collapse=",")
      m1_all<-paste0("[",m1_all,"]")
      z1<-matrix(c("samples",n1,"data",m1_all),ncol=1)

      #pie
      pie_1<-"B cells naive,B cells memory,Plasma cells,T cells CD8,T cells CD4 naive,T cells CD4 memory resting,T cells CD4 memory activated,T cells follicular helper,T cells regulatory (Tregs),T cells gamma delta,NK cells resting,NK cells activated,Monocytes,Macrophages M0,Macrophages M1,Macrophages M2,Dendritic cells resting,Dendritic cells activated,Mast cells resting,Mast cells activated,Eosinophils,Neutrophils"
      pie_1<-paste0("[",pie_1,"]")
      pie_3<-paste0(color_lm22,collapse = ",")
      for(i in 1:dim(Results)[1]){
        pie_2<-paste0(Results[i,1:22],collapse = ",")
        pie_2<-paste0("[",pie_2,"]")

        z1<-matrix(c(pie_1,pie_2,pie_3),ncol=1)

      }

    }else if(sample=="single"){

      #radar one file
      m<-"Sample"
      for(j in 1:22)
      {
        m<- paste(m,Results_radar[j],sep = ",")
      }
      z<-matrix(c(n,m),ncol=1)

      #stackplot
      n1<-paste0('["',rownames(Results),'"]')
      m1_all<-list()
      for(i in 1:22){
        m1<-paste0('{"name":"',colnames(Results)[i],'","color":"',color_lm22[i],'","data":[',Results[1,i],"]}")
        m1_all[i]<-m1
      }
      #stackplot outpur
      m1_all<-paste0(unlist(m1_all),collapse=",")
      m1_all<-paste0("[",m1_all,"]")
      z1<-matrix(c("samples",n1,"data",m1_all),ncol=1)


      #pie
      pie_1<-"B cells naive,B cells memory,Plasma cells,T cells CD8,T cells CD4 naive,T cells CD4 memory resting,T cells CD4 memory activated,T cells follicular helper,T cells regulatory (Tregs),T cells gamma delta,NK cells resting,NK cells activated,Monocytes,Macrophages M0,Macrophages M1,Macrophages M2,Dendritic cells resting,Dendritic cells activated,Mast cells resting,Mast cells activated,Eosinophils,Neutrophils"
      pie_1<-paste0("[",pie_1,"]")
      pie_3<-paste0(color_lm22,collapse = ",")
      pie_2<-paste0(Results[1:22],collapse = ",")
      pie_2<-paste0("[",pie_2,"]")

      z1<-matrix(c(pie_1,pie_2,pie_3),ncol=1)

    }

  }else if(CHIPorRNASEQ=="RNA-seq"){

    sort_name_lm14<-c(13,3,7,4,14,9,10,6,5,1,11,12,2,8)
    Results<-CIBERSORT(sig_matrix=sig_matrix,mixture_file=mixture_file,perm=perm,QN=F,sigMat="LM14",sample=sample)

    cell_type<-colnames(Results)[1:(ncol(Results)-3)]

    n1<-paste0('["',rownames(Results)[1],'"')
    color_lm14<-c("#928DBE","#77B35C","#E7D429","#F49E14","#A7D2E6","#5B9EB7","#504A8C","#E9D9C3","#EE8887","#CE482B","#E3C047","#3AB9DD","#4FAC6F","#405C70")



    Results_radar<-Results[,sort_name_lm14]
    cell_type<- cell_type[sort_name_lm14]

    n<-"name"
    for(j in 1:14){
      n<- paste(n,cell_type[j],sep = ",")
    }


    #Two processing path for parameter 'sample'
    if(sample=="multiple"){
      #radar
      for(i in 1:dim(Results_radar)[1]){
        m<-"Sample"

        for(j in 1:14)
        {
          m<- paste(m,Results_radar[i,j],sep = ",")

        }
        #radar output
        z<-matrix(c(n,m),ncol=1)

      }

      #stackplot
      for(i in 2:dim(Results)[1]){
        n1<-paste0(n1,',"',rownames(Results)[i],'"')
      }
      n1<-paste0(n1,"]")
      m1_all<-list()
      for(i in 1:14){
        m1<-paste0('{"name":"',colnames(Results)[i],'","color":"',color_lm14[i],'","data":[',Results[1,i])
        for(j in 2:dim(Results)[1]){
          m1<-paste(m1,Results[j,i],sep = ",")
        }
        m1<-paste0(m1,"]}")
        m1_all[i]<-m1
      }
      #stackplot output
      m1_all<-paste0(unlist(m1_all),collapse=",")
      m1_all<-paste0("[",m1_all,"]")
      z1<-matrix(c("samples",n1,"data",m1_all),ncol=1)



      #pie
      pie_1<-"B cells,CD4 Naive,CD4 Memory,CD8 Naive,CD8 Memory,CD8 Effector,Treg cell,Th cell,Monocytes CD16,Monocytes CD14,DC,pDC,NK,Plasma"
      pie_1<-paste0("[",pie_1,"]")
      pie_3<-paste0(color_lm14,collapse = ",")
      for(i in 1:dim(Results)[1]){
        pie_2<-paste0(Results[i,1:14],collapse = ",")
        pie_2<-paste0("[",pie_2,"]")

        z1<-matrix(c(pie_1,pie_2,pie_3),ncol=1)
      }

    }else if(sample=="single"){

      #radar
      m<-"Sample"

      for(j in 1:14)
      {
        m<- paste(m,Results_radar[j],sep = ",")

      }
      z<-matrix(c(n,m),ncol=1)

      #stackplot
      n1<-paste0('["',rownames(Results),'"]')
      m1_all<-list()
      for(i in 1:14){
        m1<-paste0('{"name":"',colnames(Results)[i],'","color":"',color_lm14[i],'","data":[',Results[1,i],"]}")
        m1_all[i]<-m1
      }
      #stackplot output
      m1_all<-paste0(unlist(m1_all),collapse=",")
      m1_all<-paste0("[",m1_all,"]")
      z1<-matrix(c("samples",n1,"data",m1_all),ncol=1)

      #pie
      pie_1<-"B cells,CD4 Naive,CD4 Memory,CD8 Naive,CD8 Memory,CD8 Effector,Treg cell,Th cell,Monocytes CD16,Monocytes CD14,DC,pDC,NK,Plasma"
      pie_1<-paste0("[",pie_1,"]")
      pie_3<-paste0(color_lm14,collapse = ",")
      pie_2<-paste0(Results[1:14],collapse = ",")
      pie_2<-paste0("[",pie_2,"]")

      z1<-matrix(c(pie_1,pie_2,pie_3),ncol=1)


    }
  }
}


#' @title CIBERSORT
#'
#' @description
#' The function is referenced by CIBERSORT_server to produce other files.
#'
#' Download this file from https://cibersort.stanford.edu/
#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Primary Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#'
#' There is a little change compared to the source code from CIBERSORT, and we already obtain permission.
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#'
#' @details
#' There are some change in parameters for our use , we need to choose diffrent signature matrix, deal with single
#' and multiple samples.
#'
#'
#' @param sig_matrix
#' Import signature matrix data,'LM22' or 'LM14'.
#'
#' @param mixture_file
#' User upload data.
#'
#' @param perm
#' Number of permutations.
#'
#' @param QN
#' Perform quantile normalization or not (TRUE/FALSE).
#'
#' @param saveDir
#' The storage path of results.
#'
#' @param sigMat
#' signature matrix name for choose,'LM22' or 'LM14'.
#'
#' @param sample
#' Choose 'multiple' or 'single' samples.
#'

CIBERSORT <- function(sig_matrix, mixture_file, perm=100, QN=TRUE,sigMat="LM22",sample="multiple"){
  # library(e1071)
  # library(parallel)
  # library(preprocessCore)

  X <- sig_matrix
  Y <- mixture_file
  X <- data.matrix(X)
  X <- X[order(rownames(X)),]

  # Two processing path for parameter 'sample'
  if(sample=="multiple"){
    Y <- data.matrix(Y)
    Y <- Y[order(rownames(Y)),]

  }else if(sample=="single"){
    Y <- as.matrix(Y)
    rownames_y<-rownames(Y)[order(rownames(Y))]
    colnames_y<-colnames(Y)
    Y <- matrix(Y[order(rownames(Y)),],ncol = 1)
    rownames(Y)<-rownames_y
    colnames(Y)<-colnames_y
  }

  P <- perm #number of permutations


  if(max(Y) < 50) {Y <- 2^Y}

  #quantile normalization of mixture file
  if(QN == TRUE){
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  if(sample=="multiple"){
    Y <- Y[YintX,]
    XintY <- Xgns %in% row.names(Y)
  }else if(sample=="single"){

    Y <- matrix(Y[YintX,],ncol=1)
    rownames(Y)<- Ygns[YintX]
    colnames(Y)<- colnames_y
    XintY <- Xgns %in% row.names(Y)

  }

  X <- X[XintY,]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if(P > 0) {nulldist <- sort(doPerm(P, X, Y)$dist)}

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")



  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- 9999

  #iterate through mixtures
  while(itor <= mixtures){

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    w <- round(result$w,3)
    mix_r <- round(result$mix_r,3)
    mix_rmse <- round(result$mix_rmse,3)

    #calculate p-value
    if(P > 0) {pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))}

    #print output
    out <- c(colnames(Y)[itor],w,round(pval,5),mix_r,mix_rmse)

    if(itor == 1) {
      output <- out
    }else {
      output <- rbind(output, out)
    }

    itor <- itor + 1

  }

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[,-1]
  obj <- obj[-1,]
  if(sample=="single"){
    obj<-matrix(obj,nrow = 1)
  }
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}


#' @title
#' Core algorithm
#'
#' @description
#' It's a subfunction of CIBERSORT.
#' SVR algorithm, the core algorithm of CIBERSORT, referenced by CIBERSORT to compute infiltration of immune cells and it can't run alone.
#'
#' @details
#' Download this file from https://cibersort.stanford.edu/
#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Primary Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' There is no change compared to the source code from CIBERSORT.
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#'
#' @param X
#' cell-specific gene expression.
#'
#' @param y
#' mixed expression per sample.
#'

CoreAlg <- function(X, y){

  #try different values of nu
  svn_itor <- 3

  res <- function(i){
    if(i==1){nus <- 0.25}
    if(i==2){nus <- 0.5}
    if(i==3){nus <- 0.75}
    model<-svm(X,y,type="nu-regression",kernel="linear",nu=nus,scale=F)
    model
  }

  if(Sys.info()['sysname'] == 'Windows') out <- mclapply(1:svn_itor, res, mc.cores=1) else
    out <- mclapply(1:svn_itor, res, mc.cores=svn_itor)

  nusvm <- rep(0,svn_itor)
  corrv <- rep(0,svn_itor)

  #do cibersort
  t <- 1
  while(t <= svn_itor) {
    weights = t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights<0)]<-0
    w<-weights/sum(weights)
    u <- sweep(X,MARGIN=2,w,'*')
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt((mean((k - y)^2)))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  #pick best model
  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]

  #get and normalize coefficients
  q <- t(model$coefs) %*% model$SV
  q[which(q<0)]<-0
  w <- (q/sum(q))

  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  newList <- list("w" = w, "mix_rmse" = mix_rmse, "mix_r" = mix_r)

}


#'
#' @title Do permutations
#'
#' @description
#' It's a subfunction of CIBERSORT.
#' To do permutations for deconvolution algorithm. Tt is the key to produce p value.
#' It referenced by CIBERSORT to compute infiltration of immune cells and it can't run alone.
#'
#' @details
#' Download this file from https://cibersort.stanford.edu/
#' CIBERSORT R script v1.03 (last updated 07-10-2015)
#' Primary Author: Aaron M. Newman, Stanford University (amnewman@stanford.edu)
#' There is no change compared to the source code from CIBERSORT.
#' License: http://cibersort.stanford.edu/CIBERSORT_License.txt
#'
#' @param perm
#' Number of permutations
#'
#' @param X
#' Cell-specific gene expression
#'
#' @param y
#' Mixed expression per sample

doPerm <- function(perm, X, Y){
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()

  while(itor <= perm){
    #print(itor)

    #random mixture
    yr <- as.numeric(Ylist[sample(length(Ylist),dim(X)[1])])

    #standardize mixture
    yr <- (yr - mean(yr)) / sd(yr)

    #run CIBERSORT core algorithm
    result <- CoreAlg(X, yr)

    mix_r <- result$mix_r

    #store correlation
    if(itor == 1) {dist <- mix_r}
    else {dist <- rbind(dist, mix_r)}

    itor <- itor + 1
  }
  newList <- list("dist" = dist)
}



