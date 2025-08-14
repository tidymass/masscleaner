# ##############svr normalization function
#' @title svrNor
#' @description svrNor
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param subject_data subject_data.
#' @param qc_data qc_data
#' @param subject_order subject_order
#' @param qc_order qc_order
#' @param multiple multiple
#' @param threads threads
#' @return result

normalize_data_svr <- function(subject_data,
                               qc_data,
                               subject_order,
                               qc_order,
                               #used data
                               multiple = 5,
                               threads = 3) {
  options(warn = -1)
  
  temp.fun <-
    function(idx,
             qc_data,
             qc_order,
             subject_data,
             subject_order,
             multiple,
             optimize_loess_span) {
      if (multiple > 1) {
        svr.reg <- e1071::svm(as.numeric(qc_data[idx, ]) ~ qc_order)
      } else {
        svr.reg <- e1071::svm(as.numeric(qc_data[idx, ]) ~ qc_order)
      }
      
      qc_data_pred <-
        summary(svr.reg)$fitted
      
      qc_nor1 <-
        as.numeric(qc_data[idx, ]) / qc_data_pred
      
      #if the predict value is 0, then set the ratio to 0
      qc_nor1[is.nan(unlist(qc_nor1))] <- 0
      qc_nor1[is.infinite(unlist(qc_nor1))] <- 0
      
      subject_data_pred <-
        predict(svr.reg, data.frame(qc_order = c(subject_order)))
      
      subject_nor1 <-
        unlist(subject_data[idx, ]) / subject_data_pred
      
      subject_nor1[is.nan(as.numeric(subject_nor1))] <-
        0
      subject_nor1[is.infinite(as.numeric(subject_nor1))] <-
        0
      subject_nor1[is.na(as.numeric(subject_nor1))] <-
        0
      subject_nor1[which(as.numeric(subject_nor1) < 0)] <-
        0
      names(qc_nor1) <- colnames(qc_data)
      return_result <- list(qc_nor1, subject_nor1)
      return(return_result)
    }
  
  
  peak_index <- seq_len(nrow(qc_data))
  
  if (masstools::get_os() == "windows") {
    bpparam <-
      BiocParallel::SnowParam(workers = threads,
                              progressbar = TRUE)
  } else{
    bpparam <-
      BiocParallel::MulticoreParam(workers = threads,
                                   progressbar = TRUE)
  }
  
  data_nor <-
    BiocParallel::bplapply(
      peak_index,
      FUN = temp.fun,
      BPPARAM = bpparam,
      qc_data = qc_data,
      qc_order = qc_order,
      subject_data = subject_data,
      subject_order = subject_order,
      multiple = multiple
    )
  
  qc_data_nor <-
    data_nor %>%
    purrr::map(function(x) {
      x[[1]]
    }) %>%
    dplyr::bind_rows()
  
  subject_data_nor <-
    data_nor %>%
    purrr::map(function(x) {
      x[[2]]
    }) %>%
    dplyr::bind_rows()
  
  qc_median <- apply(qc_data, 1, median)
  
  qc_data_nor <- qc_median * qc_data_nor
  subject_data_nor <- qc_median * subject_data_nor
  rownames(qc_data_nor) <- rownames(qc_data)
  rownames(subject_data_nor) <- rownames(qc_data)
  return_result <- list(qc_data_nor, subject_data_nor)
  message(crayon::green("SVR normalization is done."))
  return(return_result)
}
