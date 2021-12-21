#' @title integrate_data
#' @description Integrate different batch datasets together.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object object.
#' @param method Which method do you want to use?
#' subject_mean or qc_mean, default is qc_mean.
#' @return Return a mass_dataset which has been integrated.
#' @details Data integration is a necessary step for multiple batch dataset.
#' \href{https://www.readcube.com/library/fe13374b-5bc9-4c61-9b7f-6a354690947e:abe41368-d08d-4806-871f-3aa035d21743}{Dunn's}
#' method has been used in this function.
#' @export
#' @examples 
#' data("object1", package = "demodata")
#' data("object2", package = "demodata")
#' 
#' object =
#' align_batch(x = object1, y = object2, return_index = FALSE)
#' 
#' object = impute_mv(object = object, method = "zero")
#' 
#' get_mv_number(object)
#' 
#' new_object =
#' integrate_data(object = object, method = "qc_mean")
#' 
#' massdataset::intensity_plot(object = object,
#'                             variable_index = 4,
#'                             order_by = "injection.order",color_by = "batch")
#' 
#' massdataset::intensity_plot(object = new_object,
#'                             variable_index = 2,
#'                             order_by = "injection.order",color_by = "batch")

integrate_data = function(object,
                          method = c("qc_mean",
                                     "qc_median",
                                     "subject_mean",
                                     "subject_median")) {
  method <- match.arg(method)
  
  massdataset::check_object_class(object = object, class = "mass_dataset")
  
  check_result =
    check_for_data_integration(object = object, method = method)
  
  if (check_result == "warning: only one batch.") {
    warning(check_result)
    return(object)
  }
  
  if (check_result == "warning: No batch information in object.") {
    warning(check_result)
    return(object)
  }
  
  if (length(grep("error", check_result)) > 1) {
    stop(check_result)
  }
  
  if (sum(is.na(object@expression_data)) > 0) {
    stop("Please impute MV first.\n")
  }
  
  expression_data =
    object@expression_data
  
  sample_info =
    object@sample_info
  
  split_expression_data =
    purrr::map(
      unique(sample_info$batch),
      .f = function(temp_batch) {
        temp_sample_info =
          sample_info %>%
          dplyr::filter(batch == temp_batch)
        expression_data[, temp_sample_info$sample_id]
      }
    )
  
  ####correct subject samples
  correct_factor =
    purrr::map(split_expression_data, function(x) {
      qc_idx =
        which(colnames(x) %in% sample_info$sample_id[sample_info$class == "QC"])
      subject_idx =
        which(colnames(x) %in% sample_info$sample_id[sample_info$class == "Subject"])
      
      if (method == "qc_mean") {
        return(apply(x[, qc_idx], 1, mean))
      }
      
      if (method == "qc_median") {
        return(apply(x[, qc_idx], 1, median))
      }
      
      if (method == "subject_mean") {
        return(apply(x[, subject_idx], 1, mean))
      }
      
      if (method == "subject_median") {
        return(apply(x[, subject_idx], 1, median))
      }
      
    })
  
  correct_factor <-
    lapply(correct_factor, function(x) {
      correct_factor[[1]] / x
    })
  
  correct_factor <-
    correct_factor %>%
    lapply(function(x) {
      x[is.na(x)] <- 1
      x[is.infinite(x)] <- 1
      x
    })
  
  split_expression_data =
    purrr::map2(split_expression_data, correct_factor, function(x, y) {
      x * y
    })
  
  new_expression_data <-
    split_expression_data %>%
    dplyr::bind_cols()
  
  new_expression_data =
    new_expression_data[, object@sample_info$sample_id]
  
  object@expression_data = new_expression_data
  
  process_info =
    object@process_info
  
  parameter <- new(
    Class = "tidymass_parameter",
    pacakge_name = "masscleaner",
    function_name = "integrate_data()",
    parameter = list(method = method),
    time = Sys.time()
  )
  
  ###process_info
  if (all(names(process_info) != "integrate_data")) {
    process_info$integrate_data = parameter
  } else{
    process_info$integrate_data = c(process_info$integrate_data,
                                    parameter)
  }
  
  object@process_info = process_info
  
  return(object)
  
}
