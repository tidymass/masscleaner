#' @title Normalize data using total.
#' @description Normalize data using total.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param x expression data, row is variables, column is sample.
#' @param keep_scale keep scale or not.
#' @return Normalized expression data.

normalize_data_total <- function(x, keep_scale = TRUE) {
  new_x =
    x %>%
    apply(2, function(x) {
      x / sum(x)
    }) %>%
    as.data.frame()
  
  if (keep_scale) {
    max_total = max(apply(x, 2, sum))
    new_x = new_x * max_total
  }
  return(new_x)
}


#' @title Normalize data using total.
#' @description Normalize data using total.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param x expression data, row is variables, column is sample.
#' @param keep_scale keep scale or not.
#' @return Normalized expression data.

normalize_data_median <- function(x, keep_scale = TRUE) {
  new_x =
    x %>%
    apply(2, function(x) {
      x / median(x)
    }) %>%
    as.data.frame()
  
  if (keep_scale) {
    max_median = max(apply(x, 2, median))
    new_x = new_x * max_median
  }
  return(new_x)
}



#' @title Normalize data using total.
#' @description Normalize data using total.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param x expression data, row is variables, column is sample.
#' @param keep_scale keep scale or not.
#' @return Normalized expression data.

normalize_data_mean <- function(x, keep_scale = TRUE) {
  new_x =
    x %>%
    apply(2, function(x) {
      x / mean(x)
    }) %>%
    as.data.frame()
  
  if (keep_scale) {
    max_mean = max(apply(x, 2, mean))
    new_x = new_x * max_mean
  }
  return(new_x)
}