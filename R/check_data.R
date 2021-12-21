#' @title check_for_qc_normalization
#' @description check_for_qc_normalization
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object mass_dataset
#' @return error information
#' @export

check_for_qc_normalization = function(object) {
  ###QC samples
  if (all(unique(object@sample_info$class) != "QC")) {
    return("error: No QC samples in object, please check.")
  }
  
  ##injection.order
  if (all(colnames(object@sample_info) != "injection.order")) {
    return("error: No injection.order samples in object, please check.")
  }
  
  ###has batch or not
  if (all(colnames(object@sample_info) != "batch")) {
    object@sample_info$injection.order = 1
  }
  
  ###check injection.order information
  sample_info =
    object@sample_info
  
  purrr::map(
    unique(sample_info$batch),
    .f = function(batch_idx) {
      temp_sample_info =
        sample_info %>%
        dplyr::filter(batch == batch_idx)
      
      qc_order = temp_sample_info %>%
        dplyr::filter(class == "QC") %>%
        dplyr::pull(injection.order)
      
      if (min(qc_order) != 1) {
        return(paste("error: batch", batch_idx, "QC injection order is not from 1"))
      }
      
      if (max(qc_order) != max(sample_info$injection.order)) {
        return(paste(
          "error: batch",
          batch_idx,
          "QC injection order is not end last one"
        ))
      }
    }
  )
  return("ok")
}





#' @title check_for_data_integration
#' @description check_for_data_integration
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object mass_dataset
#' @param method method for integrate_data
#' @return error information
#' @export

check_for_data_integration = function(object,
                                      method = c("qc_mean",
                                                 "qc_median",
                                                 "subject_mean",
                                                 "subject_median")) {
  method <- match.arg(method)
  
  ###batch
  if (all(colnames(object@sample_info) != "batch")) {
    return("warning: No batch information in object.")
  } else{
    if (length(unique(object@sample_info$batch)) == 1) {
      return("warning: only one batch.")
    }
  }
  
  ###QC samples
  if (method == "qc_mean" |
      method == "qc_median") {
    if (all(unique(object@sample_info$class) != "QC")) {
      return("error: No QC samples in object.")
    }
    
    purrr::map(
      unique(object@sample_info$batch),
      .f = function(temp_batch) {
        number =
          object@sample_info %>%
          dplyr::filter(batch == temp_batch & class == "QC") %>%
          nrow()
        if (number == 0) {
          return(paste("error: batch", temp_batch, "has no QC."))
        }
      }
    )
  }
  
  
  ###Subject samples
  if (method == "subject_mean" |
      method == "subject_median") {
    if (all(unique(object@sample_info$class) != "Subject")) {
      return("error: No Subject samples in object.")
    }
    
    purrr::map(
      unique(object@sample_info$batch),
      .f = function(temp_batch) {
        number =
          object@sample_info %>%
          dplyr::filter(batch == temp_batch & class == "Subject") %>%
          nrow()
        if (number == 0) {
          return(paste("error: batch", temp_batch, "has no Subject"))
        }
      }
    )
  }
  
  return("ok")
  
}