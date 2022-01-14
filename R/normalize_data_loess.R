####LOESS normalization function
#' @title normalize_data_loess
#' @description normalize_data_loess
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param subject_data subject_data row is variable, column is sample
#' @param qc_data subject_data row is variable, column is sample
#' @param subject_order subject_order
#' @param qc_order qc_order
#' @param optimization optimization
#' @param begin begin see ?loess
#' @param end end ?loess
#' @param step step ?loess
#' @param threads threads
#' @return result
normalize_data_loess <- function(subject_data,
                                 qc_data,
                                 subject_order,
                                 qc_order,
                                 optimization = TRUE,
                                 begin = 0.5,
                                 end = 1,
                                 step = 0.2,
                                 threads = 4) {
  message(crayon::green("LOESS normalization...\n"))
  
  temp.fun <-
    function(idx,
             qc_data,
             qc_order,
             subject_data,
             subject_order,
             optimization = TRUE,
             begin,
             end,
             step,
             optimize_loess_span) {
      if (optimization) {
        para <- optimize_loess_span(
          x = qc_order,
          y = as.numeric(qc_data[idx,]),
          degree_range = c(1, 2),
          span_range = seq(begin, end, step)
        ) %>%
          dplyr::filter(rmse == min(rmse)) %>%
          head()
        
        loess.reg <-
          loess(
            formula = as.numeric(qc_data[idx, ]) ~ qc_order,
            span = para$span[1],
            degree = para$degree[1]
          )
      } else {
        loess.reg <- loess(as.numeric(qc_data[idx, ]) ~ qc_order)
      }
      
      qc_data_pred <-
        summary(loess.reg)$fitted
      
      qc_nor1 <-
        as.numeric(qc_data[idx, ]) / qc_data_pred
      
      #if the predict value is 0, then set the ratio to 0
      qc_nor1[is.nan(unlist(qc_nor1))] <- 0
      qc_nor1[is.infinite(unlist(qc_nor1))] <- 0
      
      subject_data_pred <-
        predict(loess.reg, data.frame(qc_order = c(subject_order)))
      
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
  
  if (tinytools::get_os() == "windows") {
    bpparam <-
      BiocParallel::SnowParam(workers = threads,
                              progressbar = TRUE)
  } else{
    bpparam <- BiocParallel::MulticoreParam(workers = threads,
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
      optimization = optimization,
      begin = begin,
      end = end,
      step = step,
      optimize_loess_span = optimize_loess_span
    )

  qc_data_nor <-
    data_nor %>% 
    purrr::map(function(x){
      x[[1]]
    }) %>% 
    dplyr::bind_rows()
  
  subject_data_nor <-
    data_nor %>% 
    purrr::map(function(x){
      x[[2]]
    }) %>% 
    dplyr::bind_rows()

  qc_median <- apply(qc_data, 1, median)
  
  qc_data_nor <- qc_median * qc_data_nor
  subject_data_nor <- qc_median * subject_data_nor
  rownames(qc_data_nor) <- rownames(qc_data)
  rownames(subject_data_nor) <- rownames(qc_data)
  return_result <- list(qc_data_nor, subject_data_nor)
  message(crayon::green("LOESS normalization is done.\n"))
  return(return_result)
}


#' @title optimize_loess_span
#' @description optimize_loess_span
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param x x
#' @param y y
#' @param degree_range numeric vector
#' @param span_range numeric vector
#' @return optimization result
#' @export
optimize_loess_span <-
  function(x,
           y,
           degree_range = c(1, 2),
           span_range = seq(0.2, 0.6, 0.1)) {
    span_rmse <-
      purrr::map(degree_range, function(degree) {
        purrr::map(span_range, function(span) {
          temp_data <-
            data.frame(x, y)
          
          prediction <-
            purrr::map(
              2:(nrow(temp_data) - 1),
              .f = function(idx) {
                temp_result <-
                  loess(
                    formula = y ~ x,
                    data = temp_data[-idx,],
                    span = span,
                    degree = degree
                  )
                prediction <-
                  try(predict(object = temp_result,
                              newdata = temp_data[idx, -2, drop = FALSE]))
                
                # if (class(prediction)[1] == "try-error") {
                if (is(prediction, class2 = "try-error")) {
                  data.frame(real = temp_data$y[idx],
                             prediction = NA)
                } else{
                  data.frame(real = temp_data$y[idx],
                             prediction = as.numeric(prediction))
                }
              }
            ) %>%
            dplyr::bind_rows()
          
          if (all(is.na(prediction$prediction))) {
            temp_rmse <- NA
          } else{
            temp_rmse <- sqrt(sum((
              prediction$real - prediction$prediction
            ) ^ 2) / nrow(prediction))
          }
          
          data.frame(span = span,
                     degree = degree,
                     rmse = temp_rmse)
        }) %>%
          dplyr::bind_rows()
      }) %>%
      dplyr::bind_rows()
    
    # plot =
    #   data.frame(x, y) %>%
    #   ggplot(aes(x, y)) +
    #   geom_point(size = 5)
    
    span_rmse <-
      span_rmse %>%
      dplyr::filter(!is.na(rmse))
    #
    # idx = which.min(span_rmse$rmse)
    # # for(i in seq_len(nrow(span_rmse))){
    # plot =
    #   plot +
    #   geom_smooth(
    #     se = FALSE,
    #     span = span_rmse$span[idx],
    #     color = ggsci::pal_lancet()(n = 9)[idx]
    #   )
    # # }
    #
    # plot =
    #   plot +
    #   ggplot2::ggtitle(label = paste("Span: ", span_rmse$span[idx])) +
    #   theme(title = element_text(colour = ggsci::pal_lancet()(n = 9)[idx]))
    
    return(span_rmse)
  }