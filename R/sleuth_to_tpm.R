#' @name sleuth_to_tpm
#'
#' @title Convert a sleuth object to TPM
#'
#' @description Convert a sleuth object to a TPM matrix with the condition names.
#'
#' @param obj \code{sleuth} object.
#' @param log bool; should natural log be returned?
#' @param pc numeric; pseudo-counts
#'
#' @return a matrix which contains a matrix of target_ids and transcript (or gene) expression in TPMs.
#'
#' @export
#'
sleuth_to_tpm <- function(obj, log = TRUE, pc = 0.5) {
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

  # extract effective length
  eff_len <- dcast(data = obs_raw, formula = target_id ~ sample, value.var = "eff_len")
  rn <- eff_len$target_id
  eff_len[, target_id := NULL]
  eff_len <- as.matrix(eff_len)
  rownames(eff_len) <- rn
  rm(rn)

  # convert to TPM
  eff_len <- eff_len[match(rownames(est_counts), rownames(eff_len)), ]
  tpm <- apply(est_counts/eff_len, 2, function(x) 1e6 * x / sum(x))

  # aggregate to gene level
  wb <- factor(so$target_mapping[match(rownames(tpm), so$target_mapping$target_id), ]$wb)
  tpm <- apply(tpm, 2, function(x) tapply(x, wb, sum))

  # normalize
  sf <- obj$tpm_sf
  tpm <- sapply(colnames(tpm), function(i) tpm[, i] / sf[i])

  # filter
  tpm <- tpm[rownames(tpm) %in% obj$filter_df$target_id, ]
  if (log) log(tpm)
  else tpm
}
