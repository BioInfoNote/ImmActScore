



#' @@ changeStepNames
#'
#' @param name.vector vector
#' @return vector

changeStepNames <- function(name.vector){
  name.vector <- gsub(".*tep1", "Step1:Release of cancer antigens", name.vector)
  name.vector <- gsub(".*tep2", "Step2:Cancer antigen presentation", name.vector)
  name.vector <- gsub(".*tep3", "Step3:Priming and activation", name.vector)
  name.vector <- gsub(".*tep4", "Step4:Trafficking of immune cells to tumors", name.vector)
  name.vector <- gsub(".*tep5", "Step5:Infiltration of immune cells into tumors", name.vector)
  name.vector <- gsub(".*tep6", "Step6:Recognition of cancer cells by T cells", name.vector)
  name.vector <- gsub(".*tep7", "Step7:Killing of cancer cells", name.vector)
  return(name.vector)

}


#' @@ filterOutliers
#' @description Filter outliers of a matrix according to the criteria of the boxplot.
#' @param  raw.matrix A matrix consist of numeric values.
#' @returnType matrix
#' @return raw matrix without outliers

filterOutliers <- function(raw.matrix){

  print("raw matrix summary before filter outliers: ")
  print(summary(as.vector(raw.matrix)))

  out <- boxplot.stats(as.vector(raw.matrix))
  print("stats:")
  print(out$stats)

  threshold_up <- out$stats[5]
  threshold_low <- out$stats[1]

  filter.number <- length(which(raw.matrix > threshold_up))
  filter.number <- filter.number + length(which(raw.matrix < threshold_low))

  print(paste0("filter.number:  ", filter.number/length(raw.matrix)))

  raw.matrix[which(raw.matrix > threshold_up)] <- threshold_up
  raw.matrix[which(raw.matrix < threshold_low)] <- threshold_low

  return(raw.matrix)
}


