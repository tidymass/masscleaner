

#' @title new_scale
#' @description new_scale
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @param new_aes new_aes
#' @return result

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

#' @title new_scale_color
#' @description new_scale_color
#' @author Xiaotao Shen
#' \email{shenxt1990@@163.com}
#' @return result
new_scale_color <- function() {
  new_scale("colour")
}




