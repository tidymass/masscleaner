#' @title Perform Probabilistic Quotient Normalization
#' @description Perform Probabilistic Quotient Normalization
#' @param x matrix to normalize. row is variable and column is sample.
#' @param pqn_reference normalization reference: "mean" for using the overall
#' average of variables as reference
#' or "median" (default) for using the overall median of variables as reference
#' @param pgn_reference_sample vector of index to specify
#' samples which average to use as reference
#' @return Normalized expression dataset.
#' @details First a total area normalization should be done before PQN is applied.
#' @export
#' @importFrom stats median
#' @author E. Nevedomskaya
#' @author Rico Derks
#' @references Dieterle, F., Ross, A., Schlotterbeck, G. & Senn, H. Probabilistic Quotient
#' Normalization as Robust Method to Account for Dilution of Complex Biological Mixtures.
#' Application in H1 NMR Metabonomics. Anal. Chem. 78, 4281-4290 (2006).

####This code is from Rico Derks (github @ricoderks), credit is his.

normalize_data_pqn <-
  function(x,
           pqn_reference = "median",
           pgn_reference_sample = NULL) {
    
    if (!is.null(pgn_reference_sample)) {
      if (pqn_reference == "mean") {
        mx <- as.numeric(apply(x[, pgn_reference_sample], 1 , mean))
      }
      if (pqn_reference == "median") {
        mx <- as.numeric(apply(x[, pgn_reference_sample], 1 , median))
      }
    } else {
      if (pqn_reference == "mean") {
        mx <- as.numeric(apply(x, 1 , mean))
      }
      if (pqn_reference == "median") {
        mx <- as.numeric(apply(x, 1 , median))
      }
    }
    
    # do the actual normalisation
    new_x <-
      purrr::map(
        x,
        .f = function(value) {
          as.numeric(value / median(as.numeric(value / mx)))
        }
      ) %>%
      dplyr::bind_cols() %>% 
      as.data.frame()
    
    rownames(new_x) <- rownames(x)
    return(new_x)
  }