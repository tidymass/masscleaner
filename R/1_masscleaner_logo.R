#' @title masscleaner_logo
#' @description masscleaner_logo
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @importFrom crayon yellow red green bold bgRed
#' @import ggplot2
#' @importFrom dplyr filter select pull everything distinct one_of left_join mutate bind_cols arrange
#' @importFrom tibble as_tibble enframe tibble rownames_to_column
#' @importFrom cli rule col_cyan tree
#' @importFrom utils packageVersion object.size write.csv tail head
#' @importFrom purrr map map2
#' @importFrom masstools get_os
#' @importFrom BiocParallel MulticoreParam SnowParam bplapply
#' @importFrom e1071 svm
#' @importFrom magrittr %>%
#' @importFrom impute impute.knn
#' @importFrom missForest missForest
#' @importFrom pcaMethods pca completeObs
#' @importFrom robust covRob
#' @importFrom stats coefficients lm loess median predict mad pchisq
#' @importFrom stats rgamma rt sd cor p.adjust prcomp t.test wilcox.test
#' @importFrom methods new is
#' @importFrom massdataset check_object_class get_sample_id show_sample_missing_values mutate_variable_na_freq activate_mass_dataset
#' @importFrom massdataset intensity_plot
#' @importClassesFrom massdataset mass_dataset tidymass_parameter
#' @import graphics
#' @return logo
#' @export

masscleaner_logo <- function() {
  message("Thank you for using masscleaner!")
  message("Version ", masscleaner_version, " (", update_date, ')')
  message("More information: masscleaner.tidymass.org")
  cat(
    c(
      "                           _____ _                            ",
      "                          / ____| |                           ",
      "  _ __ ___   __ _ ___ ___| |    | | ___  __ _ _ __   ___ _ __ ",
      " | '_ ` _ \\ / _` / __/ __| |    | |/ _ \\/ _` | '_ \\ / _ \\ '__|",
      " | | | | | | (_| \\__ \\__ \\ |____| |  __/ (_| | | | |  __/ |   ",
      " |_| |_| |_|\\__,_|___/___/\\_____|_|\\___|\\__,_|_| |_|\\___|_|   ",
      "                                                              ",
      "                                                              "
    ), sep = "\n")
}


masscleaner_version <- utils::packageVersion("masscleaner")
update_date = as.character(Sys.time())

#' @title get_masscleaner_version
#' @description get_masscleaner_version
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @return masscleaner_version
#' @export
#' @examples
#' get_masscleaner_version

get_masscleaner_version <- function() {
  return(masscleaner_version)
}

# library(cowsay)
# ##https://onlineasciitools.com/convert-text-to-ascii-art
# art <- readLines("logo.txt")
# dput(art)
