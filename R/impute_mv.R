#' @title Impute MV in data.
#' @description Impute MV in data.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param object A mass_dataset object.
#' @param sample_id which samples you want to impute missing value?
#' It is a index or character vector (sample_id)
#' @param method Imputation method. It
#' contains "knn", "rf" (missForest), "mean", "median", "zero", "minium",
#' "bpca" (BPCA), "svd" (SVD) and "ppca" (PPCA). Default is "knn".
#' The detial of
#' this method can be find in detail and reference paperes.
#' @param k See ?impute.knn
#' @param rowmax See ?impute.knn
#' @param colmax See ?impute.knn
#' @param maxp See ?impute.knn
#' @param rng.seed See ?impute.knn
#' @param maxiter See ?missForest
#' @param ntree See ?missForest
#' @param decreasing See ?missForest
#' @param nPcs See ?bpca
#' @param maxSteps See ?bpca
#' @param threshold See ?bpca
#' @param ... Other arguments.
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
#' object
#'
#' get_mv_number(object)
#' massdataset::get_mv_number(object, by = "sample")
#'
#' ###remove variables who have mv in more than 20% QC samples
#' qc_id =
#'   object %>%
#'   activate_mass_dataset(what = "sample_info") %>%
#'   filter(class == "QC") %>%
#'   pull(sample_id)
#'
#' subject_id =
#'   object %>%
#'   activate_mass_dataset(what = "sample_info") %>%
#'   filter(class == "Subject") %>%
#'   pull(sample_id)
#'
#' object =
#'   object %>%
#'   mutate_variable_na_freq(according_to_samples = qc_id) %>%
#'   mutate_variable_na_freq(according_to_samples = subject_id) %>%
#'   activate_mass_dataset(what = "variable_info") %>%
#'   filter(na_freq < 0.2 & na_freq.1 < 0.5)
#'
#' ###remove samples with MV > 50% except Blank samples
#' object =
#'   filter_samples(
#'     object = object,
#'     flist = function(x) {
#'       sum(is.na(x)) / nrow(object) < 0.5
#'     },
#'     apply_to = c(qc_id, subject_id),
#'     prune = TRUE
#'   )
#'
#' blank_id =
#'   object %>%
#'   activate_mass_dataset(what = "sample_info") %>%
#'   filter(class == "Blank") %>%
#'   pull(sample_id)
#'
#' object1 =
#'   impute_mv(object = object,
#'             sample_id = blank_id,
#'             method = "zero")
#'
#' object1 %>%
#'   activate_mass_dataset(what = "expression_data") %>%
#'   select(dplyr::contains("Blank")) %>%
#'   extract_expression_data() %>%
#'   head()
#'
#' object2 =
#'   impute_mv(object = object,
#'             sample_id = subject_id,
#'             method = "knn")
#'
#' object2 %>%
#'   activate_mass_dataset(what = "sample_info") %>%
#'   filter(class == "Subject") %>%
#'   extract_expression_data() %>%
#'   head()

####impute missing values
impute_mv <-
  function(object,
           sample_id,
           method = c("knn",
                      "rf",
                      "mean",
                      "median",
                      "zero",
                      "minimum",
                      "bpca",
                      "svdImpute",
                      "ppca"),
           k = 10,
           rowmax = 0.5,
           colmax = 0.8,
           maxp = 1500,
           rng.seed = 362436069,
           # missForest parameters
           maxiter = 10,
           ntree = 100,
           decreasing = FALSE,
           #BPCA PPCA, and SVD parameters
           nPcs = 2,
           maxSteps = 100,
           threshold = 1e-04,
           ...) {
    options(warn = -1)
    massdataset::check_object_class(object = object, class = "mass_dataset")
    
    method <- match.arg(method)
    
    if (missing(sample_id)) {
      sample_id <- massdataset::get_sample_id(object = object)
    } else{
      if (!any(sample_id %in% massdataset::get_sample_id(object = object))) {
        stop("some sample_ids are not in expression_data, please check.\n")
      }
    }
    
    #### MV imputation
    expression_data <-
      massdataset::extract_expression_data(object)
    
    expression_data1 <-
      expression_data[, sample_id]
    
    if (sum(is.na(expression_data1)) == 0) {
      return(object)
    }
    
    ## KNN method
    if (method == "knn") {
      if (exists(".Random.seed")) {
        rm(.Random.seed)
      }
      
      expression_data1 <-
        impute::impute.knn(
          data = as.matrix(expression_data1),
          k = k,
          rowmax = rowmax,
          colmax = colmax,
          maxp = maxp,
          rng.seed = rng.seed
        )
      
      expression_data1 <-
        expression_data1$data
      expression_data[, sample_id] <-
        expression_data1
      object@expression_data <-
        expression_data
      
      ##processing information
      ####add parameters
      process_info <- object@process_info
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "impute_mv()",
        parameter = list(
          method = method,
          rowmax = rowmax,
          colmax = colmax,
          maxp = maxp,
          rng.seed = rng.seed,
          sample_id = sample_id
        ),
        time = Sys.time()
      )
      
      if (all(names(process_info) != "impute_mv")) {
        process_info$impute_mv <- parameter
      } else{
        process_info$impute_mv <- c(process_info$impute_mv,
                                    parameter)
      }
      
      object@process_info <- process_info
      return(object)
    }
    
    ## rf method
    if (method == "rf") {
      expression_data1 <-
        missForest::missForest(
          t(expression_data1),
          maxiter = maxiter,
          ntree = ntree,
          decreasing = decreasing,
          ...
        )
      expression_data1 <- as.data.frame(t(expression_data1$ximp))
      
      expression_data[, sample_id] <-
        expression_data1
      
      object@expression_data <-
        expression_data
      
      ##processing information
      ####add parameters
      process_info <- object@process_info
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "impute_mv()",
        parameter = list(
          method = method,
          maxiter = maxiter,
          ntree = ntree,
          decreasing = decreasing,
          sample_id = sample_id
        ),
        time = Sys.time()
      )
      
      if (all(names(process_info) != "impute_mv")) {
        process_info$impute_mv <- parameter
      } else{
        process_info$impute_mv <- c(process_info$impute_mv,
                                    parameter)
      }
      
      object@process_info <- process_info
      return(object)
    }
    
    ## mean imputation
    if (method %in% c("mean", "median", "zero", "minimum")) {
      expression_data1 <-
        apply(expression_data1, 1, function(x) {
          x <- as.numeric(x)
          if (all(!is.na(x))) {
            return(x)
          } else{
            x[is.na(x)] <-
              switch(
                EXPR = method,
                mean = mean(x, na.rm = TRUE),
                median = median(x, na.rm = TRUE),
                zero = 0,
                minimum = min(x, na.rm = TRUE)
              )
            x
          }
        }) %>%
        t() %>%
        as.data.frame()
      colnames(expression_data1) <- sample_id
      
      expression_data[, sample_id] <-
        expression_data1
      
      object@expression_data <-
        expression_data
      
      ##processing information
      ####add parameters
      process_info <- object@process_info
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "impute_mv()",
        parameter = list(method = method),
        time = Sys.time()
      )
      
      if (all(names(process_info) != "impute_mv")) {
        process_info$impute_mv <- parameter
      } else{
        process_info$impute_mv <- c(process_info$impute_mv,
                                    parameter)
      }
      
      object@process_info <- process_info
      return(object)
    }
    
    ##BPCA, SVD and ppca
    if (method %in% c("bpca", "svdImpute", "ppca")) {
      expression_data1 <- pcaMethods::pca(
        t(expression_data1),
        method = method,
        nPcs = nPcs,
        maxSteps = maxSteps,
        threshold = threshold
      )
      
      expression_data1 <-
        as.data.frame(t(pcaMethods::completeObs(expression_data1)))
      
      expression_data[, sample_id] <-
        expression_data1
      
      object@expression_data <-
        expression_data
      
      ##processing information
      ####add parameters
      process_info <- object@process_info
      
      parameter <- new(
        Class = "tidymass_parameter",
        pacakge_name = "masscleaner",
        function_name = "impute_mv()",
        parameter = list(
          method = method,
          nPcs = nPcs,
          maxSteps = maxSteps,
          threshold = threshold
        ),
        time = Sys.time()
      )
      
      if (all(names(process_info) != "impute_mv")) {
        process_info$impute_mv <- parameter
      } else{
        process_info$impute_mv <- c(process_info$impute_mv,
                                    parameter)
      }
      
      object@process_info <- process_info
      return(object)
    }
  }
