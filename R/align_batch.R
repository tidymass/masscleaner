#' @title align_batch
#' @description Align different batch peaks tables.
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param x A mass_dataset object
#' @param y A mass_dataset object
#' @param combine.mz.tol m/z tolerance for batch alignment,
#' default is 25 ppm.
#' @param combine.rt.tol RT tolerance for batch alignment,
#' default is 30 seconds.
#' @param use.int.tol Whether use intensity match for batch alignment.
#' @param return_index return index or new object.
#' @return A index table or a new mass_dataset object.
#' @export
#' @examples
#'\dontrun{
#' data(object1, package = "demodata")
#' data(object2, package = "demodata")
#'
#' object1
#' object2
#'
#' x = object1
#' y = object2
#'
#' match_result =
#'   align_batch(x = object1, y = object2, return_index = TRUE)
#'
#' head(match_result)
#'
#' new_object =
#'   align_batch(x = object1, y = object2, return_index = FALSE)
#'
#' new_object
#' }

align_batch <- function(x,
                        y,
                        combine.mz.tol = 25,
                        combine.rt.tol = 30,
                        use.int.tol = FALSE,
                        return_index = FALSE) {
  massdataset::check_object_class(object = x, class = "mass_dataset")
  massdataset::check_object_class(object = y, class = "mass_dataset")
  
  message("Rough aligning...")
  
  rough_match_result <- rough_align(
    peak.table = list(
      cbind(x@variable_info[, c("variable_id", "mz", "rt")],
            x@expression_data),
      cbind(y@variable_info[, c("variable_id", "mz", "rt")],
            y@expression_data)
    ),
    combine.mz.tol = combine.mz.tol,
    combine.rt.tol = combine.rt.tol
  )
  
  message("Accurate aligning...")
  
  accurate_match_result <-
    accurate_align(
      peak.table = list(
        cbind(x@variable_info[, c("variable_id", "mz", "rt")],
              x@expression_data),
        cbind(y@variable_info[, c("variable_id", "mz", "rt")],
              y@expression_data)
      ),
      simple.data = rough_match_result,
      use.int.tol = use.int.tol
    )
  if (return_index) {
    return(accurate_match_result)
  } else{
    x <- x[accurate_match_result$Index1,]
    y <- y[accurate_match_result$Index2,]
    
    y@variable_info$variable_id <- x@variable_info$variable_id
    rownames(y@expression_data) <- x@variable_info$variable_id
    if (length(intersect(colnames(x), colnames(y))) > 0) {
      warning("Overlap sample IDs in x and y.\n")
    }
    return(cbind(x, y))
  }
}


#' @title rough_align
#' @description rough_align
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param peak.table peak.table
#' @param combine.mz.tol combine.mz.tol
#' @param combine.rt.tol combine.rt.tol
#' @return result

rough_align <- function(peak.table,
                        combine.mz.tol = 25,
                        combine.rt.tol = 30) {
  if (length(peak.table) == 1) {
    return(peak.table[[1]])
  }
  
  batch1 <- peak.table[[1]]
  batch2 <- peak.table[[2]]
  
  ###generate sinplifued datasets
  simple.batch1 <- simply_data(
    data = batch1,
    combine.mz.tol = combine.mz.tol,
    combine.rt.tol = combine.rt.tol
  )
  
  simple.batch2 <- simply_data(
    data = batch2,
    combine.mz.tol = combine.mz.tol,
    combine.rt.tol = combine.rt.tol
  )
  
  rm(list = c("batch1", "batch2"))
  
  ### rough matching
  data1 <- simple.batch1[, c(2:3)]
  data2 <- simple.batch2[, c(2:3)]
  
  match.result <- masstools::mz_rt_match(
    data1 = as.matrix(data1),
    data2 = as.matrix(data2),
    mz.tol = combine.mz.tol,
    rt.tol = combine.rt.tol,
    rt.error.type = "abs"
  )
  
  rm(list = c("data1", "data2"))
  
  simple.batch1 <- simple.batch1[match.result[, 1],]
  simple.batch2 <- simple.batch2[match.result[, 2],]
  
  rm("match.result")
  
  
  simple.data <- list(simple.batch1, simple.batch2)
  return(simple.data)
}


#' @title simply_data
#' @description simply_data
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data data
#' @param combine.mz.tol combine.mz.tol
#' @param combine.rt.tol combine.rt.tol
#' @return result

simply_data <- function(data,
                        combine.mz.tol = 5,
                        combine.rt.tol = 30) {
  data <- data[order(data$mz),]
  variable_id <- data$variable_id
  mz <- data$mz
  rt <- data$rt
  int <- apply(data[, -c(seq_len(3))], 1, function(x)
    mean(x, na.rm = TRUE))
  names(mz) <- names(rt) <- names(int) <- variable_id
  
  #group peaks according to mz and RT
  
  all.index <- seq_len(length(variable_id))
  remain.idx <- vector(mode = "list", length = length(variable_id))
  for (idx in seq_len(length(variable_id))) {
    if (all(idx != all.index))
      next()
    temp.mz <- mz[idx]
    temp.rt <- rt[idx]
    mz.left <- temp.mz - combine.mz.tol * temp.mz / 10 ^ 6
    mz.right <- temp.mz + combine.mz.tol * temp.mz / 10 ^ 6
    rt.left <- temp.rt - combine.rt.tol
    rt.right <- temp.rt + combine.rt.tol
    
    temp.idx <- which(mz > mz.left & mz <=  mz.right &
                        rt > rt.left & rt <= rt.right)
    if (length(temp.idx) == 1) {
      remain.idx[[idx]] <- temp.idx
    } else{
      max.idx <- temp.idx[which.max(int[temp.idx])]
      remain.idx[[idx]] <- max.idx
      all.index <- setdiff(all.index, setdiff(temp.idx, max.idx))
    }
    
  }
  
  
  remain.idx <- unique(unlist(remain.idx))
  
  simple.data <- data[remain.idx,]
  rm(list = c("data", "mz", "rt", "variable_id", "int"))
  return(simple.data)
}


#' @title baMZplot
#' @description baMZplot
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param simple.data data
#' @return result

baMZplot <- function(simple.data) {
  my.theme <- ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 18),
      axis.title.y = ggplot2::element_text(size = 18)
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 15),
      axis.text.y = ggplot2::element_text(size = 15)
    ) +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) +
    ggplot2::theme(legend.text = ggplot2::element_text(size = 10))
  
  mz.error <-
    (simple.data[[2]]$mz - simple.data[[1]]$mz) * 10 ^ 6 / simple.data[[2]]$mz
  mz1 <- simple.data[[1]]$mz
  mz.error.sd <- sd(abs(mz.error))
  mz.error.sd <-
    paste("m/z error standard:", round(mz.error.sd, 2), "ppm")
  temp.data <- data.frame(mz1, mz.error, stringsAsFactors = FALSE)
  mz.plot <- ggplot2::ggplot(data = temp.data,
                             ggplot2::aes(x = mz1, y = mz.error)) +
    ggplot2::geom_point() +
    my.theme +
    ggplot2::labs(x = "m/z (Batch 1)",
                  y = "m/z deviation (ppm)") +
    ggplot2::ggtitle(paste("m/z vs. m/z deviation; ", mz.error.sd)) +
    ggplot2::geom_hline(aes(yintercept = 0)) +
    ggplot2::annotate(
      geom = "text",
      x = -Inf,
      y = Inf,
      hjust = -0.2,
      vjust = 2,
      label = mz.error.sd
    )
  return(mz.plot)
}

#'
#' #' @title baRTplot
#' #' @description baRTplot
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param simple.data simple.data
#' #' @return result
#' baRTplot <- function(simple.data) {
#'   my.theme <- ggplot2::theme_bw() +
#'     ggplot2::theme(
#'       axis.title.x = ggplot2::element_text(size = 18),
#'       axis.title.y = ggplot2::element_text(size = 18)
#'     ) +
#'     ggplot2::theme(
#'       axis.text.x = ggplot2::element_text(size = 15),
#'       axis.text.y = ggplot2::element_text(size = 15)
#'     ) +
#'     ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) +
#'     ggplot2::theme(legend.text = ggplot2::element_text(size = 10))
#'
#'   rt.error <- simple.data[[2]]$rt - simple.data[[1]]$rt
#'   rt1 <- simple.data[[1]]$rt
#'   temp.data <- data.frame(rt1, rt.error, stringsAsFactors = FALSE)
#'   rt.error.sd <- sd(abs(rt.error))
#'   rt.error.sd <-
#'     paste("RT error standard:", round(rt.error.sd, 2), "second")
#'   rt.plot <- ggplot2::ggplot(data = temp.data,
#'                              ggplot2::aes(x = rt1, y = rt.error)) +
#'     ggplot2::geom_point() +
#'     my.theme +
#'     ggplot2::labs(x = "Retention time (s, Batch 1)",
#'                   y = "RT deviation (second)") +
#'     ggplot2::ggtitle(paste("RT vs. RT deviation;", rt.error.sd)) +
#'     ggplot2::geom_hline(aes(yintercept = 0)) +
#'     ggplot2::annotate(
#'       geom = "text",
#'       x = -Inf,
#'       y = Inf,
#'       hjust = -0.2,
#'       vjust = 2,
#'       label = rt.error.sd
#'     )
#'
#'   return(rt.plot)
#' }


#' #' @title baINTplot
#' #' @description baINTplot
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param simple.data simple.data
#' #' @return result
#'
#' baINTplot <- function(simple.data) {
#'   my.theme <- ggplot2::theme_bw() +
#'     ggplot2::theme(
#'       axis.title.x = ggplot2::element_text(size = 18),
#'       axis.title.y = ggplot2::element_text(size = 18)
#'     ) +
#'     ggplot2::theme(
#'       axis.text.x = ggplot2::element_text(size = 15),
#'       axis.text.y = ggplot2::element_text(size = 15)
#'     ) +
#'     ggplot2::theme(legend.title = ggplot2::element_text(size = 12)) +
#'     ggplot2::theme(legend.text = ggplot2::element_text(size = 10))
#'
#'   int1 <-
#'     log(apply(simple.data[[1]][, -c(seq_len(3))], 1, function(x)
#'       mean(x, na.rm = TRUE)) + 1, 10)
#'   int2 <-
#'     log(apply(simple.data[[2]][, -c(seq_len(3))], 1, function(x)
#'       mean(x, na.rm = TRUE)) + 1, 10)
#'   int.error <- int2 - int1
#'
#'   int.error.sd <- sd(abs(int.error))
#'   int.error.sd <-
#'     paste("Log10int error standard:", round(int.error.sd, 2))
#'   temp.data <- data.frame(int1, int.error, stringsAsFactors = FALSE)
#'   int.plot <- ggplot2::ggplot(data = temp.data,
#'                               ggplot2::aes(x = int1, y = int.error)) +
#'     ggplot2::geom_point() +
#'     my.theme +
#'     ggplot2::labs(x = "Log10intensity (mean, Batch 1)",
#'                   y = "Log10intensity deviation (second)") +
#'     ggplot2::ggtitle(paste("Log10int vs. Log10int deviation;", int.error.sd)) +
#'     ggplot2::geom_hline(aes(yintercept = 0)) +
#'     ggplot2::annotate(
#'       geom = "text",
#'       x = -Inf,
#'       y = Inf,
#'       hjust = -0.2,
#'       vjust = 2,
#'       label = int.error.sd
#'     )
#'   return(int.plot)
#' }


#' @title peak.table
#' @description peak.table
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param peak.table peak.table
#' @param simple.data simple.data
#' @param use.int.tol use.int.tol
#' @return result

accurate_align <- function(peak.table,
                           simple.data,
                           use.int.tol) {
  ##retrieve mz ,RT and int sd
  if (length(peak.table) == 1) {
    return(peak.table[[1]])
  }
  
  mz.error <-
    (simple.data[[2]]$mz - simple.data[[1]]$mz) * 10 ^ 6 / simple.data[[2]]$mz
  
  mz.error.sd <- sd(abs(mz.error))
  
  rt.error <- simple.data[[2]]$rt - simple.data[[1]]$rt
  rt.error.sd <- sd(abs(rt.error))
  
  int1 <-
    log(apply(simple.data[[1]][, -c(seq_len(3))], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  int2 <-
    log(apply(simple.data[[2]][, -c(seq_len(3))], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  int.error <- int2 - int1
  
  int.error.sd <- sd(abs(int.error))
  
  if (!use.int.tol) {
    int.error.sd <- int.error.sd * 1000000
  }
  
  ###begin alignment
  ref.batch <- peak.table[[1]]
  
  new.batch <- align_2batch(
    batch1 = peak.table[[1]],
    batch2 = peak.table[[2]],
    mz.error.sd = mz.error.sd,
    rt.error.sd = rt.error.sd,
    int.error.sd = int.error.sd,
    fold = 4,
    mz.weight = 0.4,
    rt.weight = 0.4,
    int.weight = 0.2
  )
  return(new.batch)
}


#' @title align_2batch
#' @description peak.table
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param batch1 batch1
#' @param batch2 batch2
#' @param mz.error.sd mz.error.sd
#' @param rt.error.sd rt.error.sd
#' @param int.error.sd int.error.sd
#' @param fold fold
#' @param mz.weight mz.weight
#' @param rt.weight rt.weight
#' @param int.weight int.weight
#' @return result

align_2batch <- function(batch1,
                         batch2,
                         mz.error.sd,
                         rt.error.sd,
                         int.error.sd,
                         fold = 4,
                         mz.weight = 0.4,
                         rt.weight = 0.4,
                         int.weight = 0.2) {
  data1 <- batch1[, c("mz", "rt")]
  int1 <-
    log(apply(batch1[, -c(seq_len(3))], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  data1 <- data.frame(data1, int1, stringsAsFactors = FALSE)
  colnames(data1)[3] <- "int"
  rownames(data1) <- batch1$name
  
  data2 <- batch2[, c("mz", "rt")]
  int2 <-
    log(apply(batch2[, -c(seq_len(3))], 1, function(x)
      mean(x, na.rm = TRUE)) + 1, 10)
  data2 <- data.frame(data2, int2, stringsAsFactors = FALSE)
  colnames(data2)[3] <- "int"
  rownames(data2) <- batch2$variable_id
  
  ##batch1 match batch2
  match.result1 <- MRImatch(
    data1 = data1,
    data2 = data2,
    mz.tol = mz.error.sd * fold,
    rt.tol = rt.error.sd * fold,
    rt.error.type = "abs",
    int.tol = int.error.sd * fold
  )
  
  unique.index1 <- unique(match.result1[, 1])
  
  remain.idx <- unlist(lapply(unique.index1, function(x) {
    temp.idx <- which(x == match.result1[, 1])
    if (length(temp.idx) == 1)
      return(temp.idx)
    temp.result <- match.result1[temp.idx,]
    ##score
    temp.mz.error <- temp.result[, "mz.error"]
    temp.rt.error <- temp.result[, "rt.error"]
    temp.int.error <- temp.result[, "int.error"]
    temp.mz.score <- lapply(temp.mz.error, function(x) {
      matchScore(error = x, sd = mz.error.sd)
    }) %>%
      unlist()
    
    temp.rt.score <- lapply(temp.rt.error, function(x) {
      matchScore(error = x, sd = rt.error.sd)
    }) %>%
      unlist()
    
    temp.int.score <- lapply(temp.int.error, function(x) {
      matchScore(error = x, sd = int.error.sd)
    }) %>%
      unlist()
    
    temp.score <-
      mz.weight * temp.mz.score + rt.weight * temp.rt.score +
      int.weight * temp.int.score
    return(temp.idx[which.max(temp.score)])
  }))
  
  
  match.result1 <- match.result1[remain.idx,]
  
  
  ##batch2 match batch1
  match.result2 <- MRImatch(
    data1 = data2,
    data2 = data1,
    mz.tol = mz.error.sd * fold,
    rt.tol = rt.error.sd * fold,
    rt.error.type = "abs",
    int.tol = int.error.sd * fold
  )
  unique.index1 <- unique(match.result2[, 1])
  
  remain.idx <- unlist(lapply(unique.index1, function(x) {
    temp.idx <- which(x == match.result2[, 1])
    if (length(temp.idx) == 1)
      return(temp.idx)
    temp.result <- match.result2[temp.idx,]
    ##score
    temp.mz.error <- temp.result[, "mz.error"]
    temp.rt.error <- temp.result[, "rt.error"]
    temp.int.error <- temp.result[, "int.error"]
    temp.mz.score <- lapply(temp.mz.error, function(x) {
      matchScore(error = x, sd = mz.error.sd)
    }) %>%
      unlist()
    
    temp.rt.score <- lapply(temp.rt.error, function(x) {
      matchScore(error = x, sd = rt.error.sd)
    }) %>%
      unlist()
    
    temp.int.score <- lapply(temp.int.error, function(x) {
      matchScore(error = x, sd = int.error.sd)
    }) %>%
      unlist()
    
    temp.score <-
      mz.weight * temp.mz.score + rt.weight * temp.rt.score +
      int.weight * temp.int.score
    return(temp.idx[which.max(temp.score)])
  }))
  
  
  match.result2 <- match.result2[remain.idx,]
  
  name1 <- paste(match.result1[, 1], match.result1[, 2], sep = "_")
  name2 <- paste(match.result2[, 2], match.result2[, 1], sep = "_")
  
  intersect.name <- intersect(name1, name2)
  
  temp.index <- which(name1 %in% intersect.name)
  
  match.result <- match.result1[temp.index,] %>%
    as.data.frame()
  
  match.result$variable_id1 <-
    batch1$variable_id[match.result$Index1]
  match.result$variable_id2 <-
    batch2$variable_id[match.result$Index2]
  
  match.result <-
    match.result %>%
    dplyr::select(variable_id1, variable_id2, dplyr::everything())
  
  # ##combine two batch data
  # batch1 <- batch1[match.result[, 1], ]
  # batch2 <- batch2[match.result[, 2], ]
  #
  # ##new name, mz and RT
  # new.mz <-
  #   data.frame(batch1[, 2], batch2[, 2], stringsAsFactors = FALSE)
  # new.mz <- apply(new.mz, 1, mean)
  #
  # new.rt <-
  #   data.frame(batch1[, 3], batch2[, 3], stringsAsFactors = FALSE)
  # new.rt <- apply(new.rt, 1, mean)
  #
  # new.name <-
  #   paste("M", round(new.mz), "T", round(new.rt), sep = "")
  #
  # new.name <- masstools::name_duplicated(new.name)
  #
  # return.data <- data.frame(new.name, new.mz, new.rt,
  #                           batch1[, -c(seq_len(3))],
  #                           batch2[, -c(seq_len(3))],
  #                           stringsAsFactors = FALSE)
  # rownames(return.data) <- new.name
  # rm(
  #   list = c(
  #     "batch1",
  #     "batch2",
  #     "new.mz",
  #     "new.rt",
  #     "new.name",
  #     "match.result",
  #     "match.result1",
  #     "match.result2"
  #   )
  # )
  return(match.result)
}



#' @title matchScore
#' @description matchScore
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param error error
#' @param sd sd
#' @return result
matchScore <- function(error, sd) {
  (sd / error) ^ 2
}


#' #' @title reName
#' #' @description reName
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param name name
#' #' @return result
#'
#' reName <- function(name) {
#'   temp.name <- unique(name)
#'   lapply(temp.name, function(x) {
#'     temp.idx <- which(x == name)
#'     if (length(temp.idx) > 1) {
#'       paste(name[temp.idx], seq_len(length(temp.idx)), sep = "_")
#'     }
#'   })
#' }


#' #' @title getBatchAlignmentInfo
#' #' @description getBatchAlignmentInfo
#' #' @author Xiaotao Shen
#' #' \email{shenxt1990@@outlook.com}
#' #' @param raw.data raw.data
#' #' @param rough.align.data rough.align.data
#' #' @param accurate.align.data accurate.align.data
#' #' @return result
#'
#'
#' getBatchAlignmentInfo <- function(raw.data,
#'                                   rough.align.data,
#'                                   accurate.align.data) {
#'   if (length(raw.data) == 1)
#'     return("Only one batch data.")
#'   mz.error <-
#'     (rough.align.data[[2]]$mz - rough.align.data[[1]]$mz) * 10 ^ 6 / rough.align.data[[2]]$mz
#'   mz.error.sd <- sd(abs(mz.error))
#'
#'   rt.error <- rough.align.data[[2]]$rt - rough.align.data[[1]]$rt
#'   rt.error.sd <- sd(abs(rt.error))
#'
#'   int1 <-
#'     log(apply(rough.align.data[[1]][, -c(seq_len(3))], 1, function(x)
#'       mean(x, na.rm = TRUE)) + 1, 10)
#'   int2 <-
#'     log(apply(rough.align.data[[2]][, -c(seq_len(3))], 1, function(x)
#'       mean(x, na.rm = TRUE)) + 1, 10)
#'   int.error <- int2 - int1
#'   int.error.sd <- sd(abs(int.error))
#'
#'   parameter.info <-
#'     paste(
#'       "The standard for m/z, RT and log10intensity errors are",
#'       round(mz.error.sd, 2),
#'       ",",
#'       round(rt.error.sd, 2),
#'       ",",
#'       'and',
#'       round(int.error.sd, 2),
#'       ",respectively. So the tolerances for m/z, RT and log10intensity are",
#'       4 * round(mz.error.sd, 2),
#'       ",",
#'       4 * round(rt.error.sd, 2),
#'       ",",
#'       'and',
#'       4 * round(int.error.sd, 2),
#'       ",respectively."
#'     )
#'
#'   peak.info <-
#'     paste(
#'       "There are ",
#'       length(raw.data),
#'       " batches. And the peak number are ",
#'       paste(unlist(lapply(raw.data, nrow)), collapse = ","),
#'       " respectively. ",
#'       "After batch alignment, the peak number is ",
#'       nrow(accurate.align.data),
#'       ".",
#'       sep = ""
#'     )
#'
#'   return(paste(parameter.info, peak.info))
#' }



#' @title MRImatch
#' @description MRImatch
#' @author Xiaotao Shen
#' \email{shenxt1990@@outlook.com}
#' @param data1 data1
#' @param data2 data2
#' @param mz.tol mz.tol
#' @param rt.tol rt.tol
#' @param rt.error.type rt.error.type
#' @param int.tol int.tol
#' @return result

MRImatch <- function(data1,
                     data2,
                     mz.tol,
                     #rt.tol is relative
                     rt.tol = 30,
                     rt.error.type = c("relative", "abs"),
                     int.tol = 1) {
  rt.error.type <- match.arg(rt.error.type)
  #
  if (nrow(data1) == 0 | nrow(data2) == 0) {
    result <- NULL
    return(result)
  }
  # mz1 <- as.numeric(data1[, 1])
  # rt1 <- as.numeric(data1[, 2])
  info1 <- data1[, c(1, 2, 3)]
  info1 <- apply(info1, 1, list)
  
  mz2 <- as.numeric(data2[, 1])
  rt2 <- as.numeric(data2[, 2])
  int2 <- as.numeric(data2[, 3])
  
  result <-
    info1 %>%
    purrr::map(function(x) {
      temp.mz1 <- x[[1]][[1]]
      temp.rt1 <- x[[1]][[2]]
      temp.int1 <- x[[1]][[3]]
      mz.error <- abs(temp.mz1 - mz2) * 10 ^ 6 / temp.mz1
      if (rt.error.type == "relative") {
        rt.error <- abs(temp.rt1 - rt2) * 100 / temp.rt1
      } else{
        rt.error <- abs(temp.rt1 - rt2)
      }
      
      int.error <- abs(temp.int1 - int2)
      
      j <-
        which(mz.error <= mz.tol &
                rt.error <= rt.tol & int.error <= int.tol)
      if (length(j) == 0) {
        matrix(NA, ncol = 10)
      } else{
        cbind(
          j,
          temp.mz1,
          mz2[j],
          mz.error[j],
          temp.rt1,
          rt2[j],
          rt.error[j],
          temp.int1,
          int2[j],
          int.error[j]
        )
      }
    })
  
  if (length(result) == 1) {
    result <- cbind(1, result[[1]])
  } else{
    result <- mapply(function(x, y) {
      list(cbind(x, y))
    },
    x <- seq_len(length(info1)),
    y = result)
    result <- do.call(rbind, result)
  }
  
  result <-
    matrix(result[which(!apply(result, 1, function(x)
      any(is.na(x)))),], ncol = 11)
  if (nrow(result) == 0)
    return(NULL)
  colnames(result) <-
    c(
      "Index1",
      "Index2",
      "mz1",
      "mz2",
      "mz.error",
      "rt1",
      "rt2",
      "rt.error",
      "int1",
      "int2",
      "int.error"
    )
  result <- result
}
