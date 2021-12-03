#' @title masscleaner_logo
#' @description masscleaner_logo
#' @author Xiaotao Shen
#' \email{shenxt@@stanford.edu}
#' @importFrom crayon yellow red green bold bgRed
#' @import ggplot2
#' @importFrom pbapply pblapply pboptions
#' @importFrom stringr str_split str_replace_all str_trim str_detect str_extract
#' @importFrom dplyr filter select pull everything distinct one_of left_join mutate bind_cols arrange
#' @importFrom tibble as_tibble enframe tibble rownames_to_column
#' @importFrom clisymbols symbol
#' @importFrom cli rule col_cyan tree
#' @importFrom utils packageVersion object.size write.csv tail
#' @importFrom purrr map map2
#' @importFrom plyr dlply .
#' @importFrom RColorBrewer brewer.pal
#' @importFrom readr read_csv cols
#' @importFrom readxl read_excel
#' @importFrom ggrepel geom_text_repel
#' @importFrom tinytools get_os
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @importFrom e1071 svm
#' @importFrom magrittr %>%
#' @importFrom plotly ggplotly
#' @importFrom BiocGenerics basename
#' @importFrom impute impute.knn
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca completeObs
#' @importFrom patchwork plot_layout
#' @import patchwork
#' @importFrom stats coefficients lm loess median predict 
#' @importFrom stats rgamma rt sd cor p.adjust prcomp t.test wilcox.test
#' @export

masscleaner_logo <- function(){
  cat(crayon::green("Thank you for using masscleaner!\n"))
  cat(crayon::green("Version 0.9.2 (20210312)\n"))
  cat(crayon::green("Bug fixing\n"))
  cat(crayon::green("More information can be found at https://tidymass.github.io/masscleaner/\n"))
  cat(crayon::green(
    c("                 _    __ _              ___  ", "                | |  / _| |            |__ \\ ",
      "  _ __ ___   ___| |_| |_| | _____      __ ) |", " | '_ ` _ \\ / _ \\ __|  _| |/ _ \\ \\ /\\ / // / ",
      " | | | | | |  __/ |_| | | | (_) \\ V  V // /_ ", " |_| |_| |_|\\___|\\__|_| |_|\\___/ \\_/\\_/|____|",
      "                                             ", "                                             "
    )
    
  ), sep = "\n")
}




