#' @name sleuth_to_count
#'
#' @title Convert a sleuth object to count
#'
#' @description Convert a sleuth object to a count matrix with the condition names.
#'
#' @param obj \code{sleuth} object.
#' @param log bool; should natural log be returned?
#' @param pc numeric; pseudo-counts
#'
#' @return a matrix which contains a matrix of target_ids and transcript (or gene) expression in read counts.
#'
#' @export
#'
sleuth_to_count <- function(obj, log = TRUE, pc = 0.5) {
  obs_raw <- as.data.table(obj$obs_raw)
  setkeyv(obs_raw, "target_id")

  est_counts <- dcast(data = obs_raw, formula = target_id ~ sample, value.var = "est_counts")
  rn <- est_counts$target_id
  est_counts[, target_id := NULL]
  est_counts <- as.matrix(est_counts)
  rownames(est_counts) <- rn
  rm(rn)

  # add pseudo cout
  est_counts <- est_counts + pc

  # aggregate to gene level
  wb <- factor(so$target_mapping[match(rownames(est_counts), so$target_mapping$target_id), ]$wb)
  est_counts <- apply(est_counts, 2, function(x) tapply(x, wb, sum))

  # normalize
  sf <- obj$est_counts_sf
  est_counts <- sapply(colnames(est_counts), function(i) est_counts[, i] / sf[i])

  # filter
  est_counts <- est_counts[rownames(est_counts) %in% obj$filter_df$target_id, ]
  if (log) log(est_counts)
  else est_counts
}
