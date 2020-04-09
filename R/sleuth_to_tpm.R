#' @name sleuth_to_tpm
#'
#' @title Convert a sleuth object to TPM
#'
#' @description Convert a sleuth object to a TPM matrix with the condition names.
#'
#' @inheritParams sleuth_to_count
#'
#' @return a matrix which contains a matrix of target_ids and transcript (or gene) expression in TPMs.
#'
#' @export
#'
sleuth_to_tpm <- function(obj,
                          log = TRUE,
                          normalize = TRUE,
                          pc = 0.5,
                          aggregation_column = "wb") {

  obs_raw <- as.data.table(obj$obs_raw)
  est_counts <- castToMatrix(data = obs_raw,
                      formula = target_id ~ sample,
                      value.var = "est_counts")

  # add pseudo cout
  est_counts <- est_counts + pc

  # extract effective length
  eff_len <- castToMatrix(data = obs_raw,
                          formula = target_id ~ sample,
                          value.var = "eff_len")

  # convert to TPM
  eff_len <- eff_len[match(rownames(est_counts), rownames(eff_len)), ]
  tpm <- apply(est_counts/eff_len, 2, function(x) 1e6 * x / sum(x))

  # aggregate to gene level
  ag <- factor(obj$target_mapping[match(rownames(tpm), obj$target_mapping$target_id), ][, aggregation_column])
  tpm <- apply(tpm, 2, function(x) tapply(x, ag, sum))

  # normalize
  if (normalize) {
    sf <- obj$tpm_sf
    tpm <- sapply(colnames(tpm), function(i) tpm[, i] / sf[i])
  }

  # filter
  tpm <- tpm[rownames(tpm) %in% obj$filter_df$target_id, ]
  if (log) log(tpm)
  else tpm
}
