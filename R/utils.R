
msg <- function(..., startup = FALSE) {
  if (startup) {
    if (!isTRUE(getOption("masscleaner.quiet"))) {
      packageStartupMessage(text_col(...))
    }
  } else {
    message(text_col(...))
  }
}

text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  
  theme <- rstudioapi::getThemeInfo()
  
  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)
  
}

#' List all packages in the masscleaner
#'
#' @param include_self Include masscleaner in the list?
#' @export
#' @examples
#' masscleaner_packages()
masscleaner_packages <- function(include_self = TRUE) {
  raw <- utils::packageDescription("masscleaner")$Imports
  imports <- strsplit(raw, ",")[[1]]
  parsed <- gsub("^\\s+|\\s+$", "", imports)
  names <- vapply(strsplit(parsed, "\\s+"), "[[", 1, FUN.VALUE = character(1))
  
  if (include_self) {
    names <- c(names, "masscleaner")
  }
  
  names
}

invert <- function(x) {
  if (length(x) == 0) return()
  stacked <- utils::stack(x)
  tapply(as.character(stacked$ind), stacked$values, list)
}


style_grey <- function(level, ...) {
  crayon::style(
    paste0(...),
    crayon::make_style(grDevices::grey(level), grey = TRUE)
  )
}

#' @title new_scale
#' @description new_scale
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param new_aes new_aes
#' @return result

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

#' @title new_scale_color
#' @description new_scale_color
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @return result
new_scale_color <- function() {
  new_scale("colour")
}




