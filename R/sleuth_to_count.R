#' @name sleuth_to_count
#'
#' @title Convert a sleuth object to count
#'
#' @description Convert a sleuth object to a count matrix with the condition names.
#'
#' @param obj \code{sleuth} object.
#' @param log bool; should natural log be returned?
#' @param normalize bool; should normalized counts be returned
#' @param pc numeric; pseudo-counts
#' @param aggregation_column character; aggregation column
#'
#'
#' @return a matrix which contains a matrix of target_ids and transcript (or gene) expression in read counts.
#'
#' @export
#'
sleuth_to_count <- function(obj,
                            log = TRUE,
                            normalize = TRUE,
                            pc = 0.5,
                            aggregation_column = "wb") {


  obs_raw <- as.data.table(obj$obs_raw)
  setkeyv(obs_raw, "target_id")

  est_counts <- castToMatrix(data = obs_raw,
                             formula = target_id ~ sample,
                             value.var = "est_counts")

  # add pseudo cout
  est_counts <- est_counts + pc

  # aggregate to gene level
  ag <- factor(obj$target_mapping[match(rownames(est_counts), obj$target_mapping$target_id), ][, aggregation_column])
  est_counts <- apply(est_counts, 2, function(x) tapply(x, ag, sum))

  # normalize
  if (normalize) {
    sf <- obj$est_counts_sf
    est_counts <- sapply(colnames(est_counts), function(i) est_counts[, i] / sf[i])
  }

  # filter
  est_counts <- est_counts[rownames(est_counts) %in% obj$filter_df$target_id, ]
  if (log) log(est_counts)
  else est_counts
}
