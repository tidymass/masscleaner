#' @title Detect outlier samples
#' @description Detect outlier samples. See more here:
#' \url{https://privefl.github.io/blog/detecting-outlier-samples-in-pca/}
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object A mass_dataset object.
#' @param na_percentage_cutoff na_percentage_cutoff
#' @param sd_fold_change sd_fold_change
#' @param mad_fold_change mad_fold_change
#' @param dist_p_cutoff dist_p_cutoff
#' @return A new mass_dataset object.
#' @export
#' @examples
#' library(massdataset)
#' data("expression_data")
#' data("sample_info")
#' data("variable_info")
#' object =
#'   create_mass_dataset(
#'     expression_data = expression_data,
#'     sample_info = sample_info,
#'     variable_info = variable_info
#'   )
#' 
#' object =
#'   object %>%
#'   log() %>%
#'   scale()
#' 
#' outlier_samples =
#'   object %>%
#'   detect_outlier()
#' 
#' extract_outlier_table(outlier_samples)
#' 
#' ###MV plot
#' massdataset::show_sample_missing_values(object = object,
#'                                         color_by = "class",
#'                                         percentage = TRUE)

####impute missing values
detect_outlier <- function(object,
                          na_percentage_cutoff = 0.5,
                          sd_fold_change = 6,
                          mad_fold_change = 6,
                          dist_p_cutoff = 0.05) {
  options(warn = -1)
  massdataset::check_object_class(object = object, class = "mass_dataset")
  
  ####according to NA percentage
  
  na_percentage <-
    apply(object@expression_data, 2, function(x) {
      sum(is.na(x)) / nrow(object@expression_data)
    })
  
  according_to_na <-
    na_percentage > na_percentage_cutoff
  
  according_to_na <-
    data.frame(according_to_na) %>%
    tibble::rownames_to_column(var = "sample_id")
  
  ###PCA
  if (sum(is.na(object@expression_data)) > 0) {
    warning("MVs in you object,\nwill remove variables > 50% and imputate with zero.\n")
    object <-
      object %>%
      massdataset::mutate_variable_na_freq()
    object <-
      object %>%
      massdataset::activate_mass_dataset(what = "variable_info") %>%
      dplyr::filter(na_freq < 0.5)
  }
  
  sample_info <- object@sample_info
  expression_data <- object@expression_data
  
  expression_data <-
    expression_data %>%
    apply(1, function(x) {
      x[is.na(x)] <- min(x[!is.na(x)])
      x
    }) %>%
    t()
  
  pca_object <- prcomp(x = t(as.matrix(expression_data)),
                      center = FALSE,
                      scale. = FALSE)
  
  pc <- as.data.frame(pca_object$x)
  
  # pc %>%
  #   ggplot(aes(PC1, PC2)) +
  #   geom_point()
  
  #####according to first PCs and if they out of 6 sds
  #####The standard way to detect outliers in genetics is the
  #####criterion of being “more than 6 standard deviations away from the mean”.
  according_to_pc_sd <-
    apply(pc[, c(1, 2)], 2, function(x)
      abs(x - mean(x)) > (sd_fold_change * sd(x))) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    as.data.frame() %>%
    dplyr::rename(pc1_sd = PC1,
                  pc2_sd = PC2)
  
  according_to_pc_sd$pc_sd <-
    apply(according_to_pc_sd[, -1], 1, function(x) {
      any(x)
    })
  
  according_to_pc_sd <-
    according_to_pc_sd %>%
    dplyr::select(sample_id, pc_sd)
  
  #####according to first PCs and if they out of 6 mads
  #####TNote that you might want to use median() instead of mean() and mad()
  #####instead of sd() because they are more robust estimators. This becomes
  according_to_pc_mad <-
    apply(pc[, c(1, 2)], 2, function(x)
      abs(x - median(x)) > (mad_fold_change * mad(x))) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "sample_id") %>%
    as.data.frame() %>%
    dplyr::rename(pc1_mad = PC1,
                  pc2_mad = PC2)
  
  according_to_pc_mad$pc_mad <-
    apply(according_to_pc_mad[, -1], 1, function(x) {
      any(x)
    })
  
  according_to_pc_mad <-
    according_to_pc_mad %>%
    dplyr::select(sample_id, pc_mad)
  
  #####according to distance
  ##Robust Mahalanobis distance
  ##Instead of using the infinite distance, Mahalanobis distance is
  ##a multivariate distance based on all variables (PCs here) at once.
  ##We use a robust version of this distance, which is implemented in
  ##packages {robust} and {robustbase} (Gnanadesikan and
  ##Kettenring 1972, Yohai and Zamar (1988), Maronna and Zamar (2002),
  ##Todorov, Filzmoser, and others (2009)) and that is reexported
  ## in {bigutilsr}.
  dist <- robust::covRob(data = pc[, c(seq_len(2))],
                         estim = "pairwiseGK")$dist
  
  ##This new criterion provides similar results for this data.
  ##These robust Mahalanobis distances are approximately
  ## Chi-square distributed, which enables deriving p-values of outlierness.
  
  pval <- pchisq(dist, df = 10, lower.tail = FALSE)
  # hist(pval)
  
  accordint_to_distance <-
    (pval < (dist_p_cutoff / length(dist)))  # Bonferroni correction
  
  accordint_to_distance <-
    data.frame(accordint_to_distance) %>%
    tibble::rownames_to_column(var = "sample_id")
  
  # ###Local Outlier Factor (LOF)
  # ###LOF statistic (Breunig et al. 2000) has been cited more than 4000 times.
  # ##Instead of computing a distance from the center,
  # ##it uses some local density of points. We make use of the fast
  # ##K nearest neighbours implementation of R package {nabor}
  # ##(Elseberg et al. 2012) to implement this
  # ##statistic efficiently in {bigutilsr}.
  # llof <-
  #   tryCatch(expr = bigutilsr::LOF(U = pc), error = function(e){
  #     NA
  #   })
  #
  # if(is.na(llof)[1]){
  #
  # }
  
  ####output results
  ####
  
  result <- new(
    Class = "outlier_samples",
    outlier_samples_table = list(
      according_to_na = according_to_na,
      according_to_pc_sd = according_to_pc_sd,
      according_to_pc_mad = according_to_pc_mad,
      accordint_to_distance = accordint_to_distance
    ),
    parameter = list(
      according_to_na = list(na_percentage_cutoff = na_percentage_cutoff),
      according_to_pc_sd = list(sd_fold_change = sd_fold_change),
      according_to_pc_mad = list(mad_fold_change = mad_fold_change),
      accordint_to_distance = list(dist_p_cutoff = dist_p_cutoff)
    ),
    intermediate = list(
      according_to_na = list(pc = pc),
      according_to_pc_sd = list(),
      according_to_pc_mad = list(),
      accordint_to_distance = list(dist = dist)
    ),
    ref = list(
      according_to_na = "privefl.github.io/blog/detecting-outlier-samples-in-pca",
      according_to_pc_sd = "privefl.github.io/blog/detecting-outlier-samples-in-pca",
      according_to_pc_mad = "privefl.github.io/blog/detecting-outlier-samples-in-pca",
      accordint_to_distance = "privefl.github.io/blog/detecting-outlier-samples-in-pca"
    )
  )
  
  return(result)
  
}


##S4 class for outlier samples
#' An S4 class that stores the outlier sample information.
#' @slot outlier_samples_table outlier_samples_table
#' @slot parameter parameter
#' @slot parameter parameter
#' @slot intermediate intermediate
#' @slot ref red
#' @exportClass outlier_samples
setClass(
  Class = "outlier_samples",
  representation(
    outlier_samples_table = "list",
    parameter = "list",
    intermediate = "list",
    ref = "list"
  )
)


setMethod(
  f = "show",
  signature = "outlier_samples",
  definition = function(object) {
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    cat(crayon::green("masscleaner", "\n"))
    cat(crayon::yellow(paste(rep("-", 20), collapse = ""), "\n"))
    for (i in seq_len(length(object@outlier_samples_table))) {
      cat(crayon::green(i, names(object@outlier_samples_table)[i], ": "))
      cat(sum(object@outlier_samples_table[[i]][, 2]),
          "outlier samples.\n")
      if (sum(object@outlier_samples_table[[i]][, 2]) > 0) {
        sample_id <-
          object@outlier_samples_table[[i]]$sample_id[object@outlier_samples_table[[i]][, 2]]
        if (length(sample_id) > 5) {
          sample_id <- paste(c(head(sample_id, 5), "..."), collapse = ",")
        } else{
          sample_id <- paste(head(sample_id, 5), collapse = ",")
        }
        cat(crayon::green(sample_id, ".\n"))
        
      }
    }
    cat(crayon::yellow("extract all outlier table using extract_outlier_table()\n"))
  }
)


#' @title extract_outlier_table
#' @description extract_outlier_table
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object outlier_samples class object.
#' @return A data frame.
#' @export

extract_outlier_table <- function(object) {
  object@outlier_samples_table %>%
    lapply(function(x) {
      x %>%
        tibble::column_to_rownames(var = "sample_id")
    }) %>%
    dplyr::bind_cols()
}