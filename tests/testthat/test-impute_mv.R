library(massdataset)
library(testthat)
data("expression_data")
data("sample_info")
data("variable_info")
object =
  create_mass_dataset(
    expression_data = expression_data,
    sample_info = sample_info,
    variable_info = variable_info
  )

object

sum(is.na(object))

###remove variables who have mv in more than 20% QC samples
qc_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "QC") %>%
  pull(sample_id)

subject_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "Subject") %>%
  pull(sample_id)

object =
  object %>%
  mutate_variable_na_freq(according_to_samples = qc_id) %>%
  mutate_variable_na_freq(according_to_samples = subject_id) %>%
  activate_mass_dataset(what = "variable_info") %>%
  filter(na_freq < 0.2 & na_freq.1 < 0.5)

###remove samples with MV > 50% except Blank samples
object =
  filter_samples(
    object = object,
    flist = function(x) {
      sum(is.na(x)) / nrow(object) < 0.5
    },
    apply_to = c(qc_id, subject_id),
    prune = TRUE
  )

blank_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(class == "Blank") %>%
  pull(sample_id)

test_that("impute_mv", {
  object1 <-
    impute_mv(object = object,
              sample_id = blank_id,
              method = "zero")
  
  testthat::expect_s4_class(object = object1,
                            class = "mass_dataset")
  
  object2 <-
    impute_mv(object = object,
              sample_id = subject_id,
              method = "knn")
  testthat::expect_s4_class(object = object2,
                            class = "mass_dataset")
  
  testthat::expect_error(object =  impute_mv(
    object = object@expression_data,
    sample_id = subject_id,
    method = "knn"
  ))
  
})
