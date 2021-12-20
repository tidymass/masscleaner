#' @title normalize_data
#' @description Data normalization.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object A mass_dataset object.
#' @param method Normalization method, mean, median, total svr or loess,
#' default is svr. Please see the details.
#' @param keep_scale Remain scale or not. Default is TRUE.
#' @param pqn_reference for pqn method.
#' @param begin 0.5
#' @param end 1
#' @param step 0.2
#' @param multiple multiple
#' @param threads 4
#' @export
#' @return A new mass_dataset object.
#' data("object1", package = "demodata")
#' object1 = impute_mv(object1, method = "minimum")
#' object_mean = normalize_data(object = object1, method = "mean")

normalize_data =
  function(object,
           method = c("svr", "total", "median", "mean", "pqn", "loess"),
           keep_scale = TRUE,
           pqn_reference = c("median", "mean"),
           begin = 0.5,
           end = 1,
           step = 0.2,
           multiple = 1,
           threads = 4) {
    method <- match.arg(method)
    pqn_reference = match.arg(pqn_reference)
    
    massdataset::check_object_class(object = object, class = "mass_dataset")
    
    if (method == "svr" | method == "loess") {
      check_result =
        check_for_qc_normalization(object = object)
      if (length(grep("error", check_result)) > 0) {
        stop(check_result)
      }
    }
    
    expression_data =
      object@expression_data
    
    if (sum(is.na(expression_data)) > 0) {
      stop("Please impute MV first.\n")
    }
    
    ##processing information
    ####add parameters
    process_info = object@process_info
    
    ##sample-wise methods
    if (method == "total") {
      new_expression_data =
        normalize_data_total(x = expression_data, keep_scale = keep_scale)
      
      object@expression_data = new_expression_data
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "normalize_data()",
        parameter = list(method = method,
                         keep_scale = keep_scale),
        time = Sys.time()
      )
    }
    
    if (method == "mean") {
      new_expression_data =
        normalize_data_mean(x = expression_data, keep_scale = keep_scale)
      object@expression_data = new_expression_data
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "normalize_data()",
        parameter = list(method = method,
                         keep_scale = keep_scale),
        time = Sys.time()
      )
    }
    
    if (method == "median") {
      new_expression_data =
        normalize_data_median(x = expression_data, keep_scale = keep_scale)
      object@expression_data = new_expression_data
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "normalize_data()",
        parameter = list(method = method,
                         keep_scale = keep_scale),
        time = Sys.time()
      )
    }
    
    ####-----------------------------------------------------------------------
    ##pqn (Probabilistic Quotient Normalization) method
    if (method == "pqn") {
      pgn_reference_sample =
        which(object@sample_info$class == "QC")
      if (length(pgn_reference_sample) == 0) {
        pgn_reference_sample = NULL
      }
      new_expression_data =
        normalize_data_pqn(
          x = expression_data,
          pqn_reference = pqn_reference,
          pgn_reference_sample = pgn_reference_sample
        )
      
      object@expression_data = new_expression_data
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "normalize_data()",
        parameter = list(
          method = method,
          keep_scale = keep_scale,
          pqn_reference = pqn_reference,
          pgn_reference_sample = pgn_reference_sample
        ),
        time = Sys.time()
      )
      
    }
    
    #######loess normalization
    sample_info =
      object@sample_info
    
    if (all(colnames(object@sample_info) != "batch")) {
      sample_info$batch = 1
    }
    
    if (method == "loess") {
      data_nor =
        purrr::map(unique(sample_info$batch), function(batch_idx) {
          cat(crayon::yellow("Batch", batch_idx, "...", "\n"))
          subject_id =
            sample_info %>%
            dplyr::filter(class == "Subject" &
                            batch == batch_idx) %>%
            dplyr::pull(sample_id)
          
          subject_idx = match(subject_id, sample_info$sample_id)
          
          qc_id =
            sample_info %>%
            dplyr::filter(class == "QC" & batch == batch_idx) %>%
            dplyr::pull(sample_id)
          
          qc_idx = match(qc_id, sample_info$sample_id)
          
          subject_data = expression_data[, subject_idx]
          qc_data = expression_data[, qc_idx]
          
          subject_order = sample_info$injection.order[subject_idx]
          qc_order = sample_info$injection.order[qc_idx]
          
          new_data =
            normalize_data_loess(
              subject_data = subject_data,
              qc_data = qc_data,
              subject_order = subject_order,
              qc_order = qc_order,
              optimization = optimization,
              begin = begin,
              end = end,
              step = step,
              threads = threads
            )
          
          new_data =
            new_data %>%
            dplyr::bind_cols()
          
          new_data
        })
      
      data_nor =
        data_nor %>%
        dplyr::bind_cols()
      
      new_expression_data =
        expression_data
      
      new_expression_data[, colnames(data_nor)] =
        data_nor
      
      object@expression_data = new_expression_data
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "normalize_data()",
        parameter = list(
          method = method,
          keep_scale = keep_scale,
          optimization = optimization,
          begin = begin,
          end = end,
          step = step,
          threads = threads
        ),
        time = Sys.time()
      )
      
      
    }
    
    if (method == "svr") {
      data_nor =
        purrr::map(unique(sample_info$batch), function(batch_idx) {
          cat(crayon::yellow("Batch", batch_idx, "...", "\n"))
          subject_id =
            sample_info %>%
            dplyr::filter(class == "Subject" &
                            batch == batch_idx) %>%
            dplyr::pull(sample_id)
          
          subject_idx = match(subject_id, sample_info$sample_id)
          
          qc_id =
            sample_info %>%
            dplyr::filter(class == "QC" & batch == batch_idx) %>%
            dplyr::pull(sample_id)
          
          qc_idx = match(qc_id, sample_info$sample_id)
          
          subject_data = expression_data[, subject_idx]
          qc_data = expression_data[, qc_idx]
          
          subject_order = sample_info$injection.order[subject_idx]
          qc_order = sample_info$injection.order[qc_idx]
          
          new_data =
            normalize_data_svr(
              subject_data = subject_data,
              qc_data = qc_data,
              subject_order = subject_order,
              qc_order = qc_order,
              multiple = multiple,
              threads = threads
            )
          
          new_data =
            new_data %>%
            dplyr::bind_cols()
          
          new_data
        })
      
      data_nor =
        data_nor %>%
        dplyr::bind_cols()
      
      new_expression_data =
        expression_data
      
      new_expression_data[, colnames(data_nor)] =
        data_nor
      
      object@expression_data = new_expression_data
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "normalize_data()",
        parameter = list(
          method = method,
          keep_scale = keep_scale,
          multiple = multiple,
          threads = threads
        ),
        time = Sys.time()
      )
      
    }
    
    
    ###process_info
    if (all(names(process_info) != "normalize_data")) {
      process_info$normalize_data = parameter
    } else{
      process_info$normalize_data = c(process_info$normalize_data,
                                      parameter)
    }
    
    object@process_info = process_info
    
    
    invisible(object)
  }
