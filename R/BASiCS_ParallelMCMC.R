#' @name BASiCS_ParallelMCMC
#'
#' @title Parallelized BASiCS_MCMC
#'
#' @description Parallelize BASiCS normalization implemented in \link[BASiCS]{BASiCS_MCMC}.
#'
#' @details Assumes that there are spike-ins and that these are named "ERCC".
#' Filter rule: a gene must be counted at least \code{MinCount} in \code{PropMinCount} percent of the samples.
#' For biological transcripts, this needs to be true in at least one group. For spike-ins in all the groups.
#' Uses \link[parallel]{mclapply} to parallelize.
#'
#' @param CountsList A list of unfiltered count matrices with genes as rows and samples as columns.
#' @inheritParams BASiCS::BASiCS_MCMC
#' @inheritParams BASiCS::newBASiCS_Data
#' @param MinCount integer; Minimum number of counts for a gene.
#' @param PropMinCount; numeric; Proportion of samples with minimum count.
#' @param PropMinCount; numeric; Proportion of samples with minimum count.
#' @param NCores; integer; Number of cores.
#' @param Name character; Chain identification name.
#' @param Verbose; bool; Print progress in \code{./logs} directory?
#'
#' @return A list of \link[BASiCS]{BASiCS_Chain} object.
#'
#' @seealso \link[BASiCS]{BASiCS_Chain}
#' @seealso \link[BASiCS]{BASiCS_MCMC}
#'
#' @export
#'
BASiCS_ParallelMCMC <- function(CountsList,
                                SpikeInfo,
                                N = 200000L,
                                Burn = 80000L,
                                Thin = 40L,
                                Regression = TRUE,
                                MinCount = 10L,
                                PropMinCount = 0.85,
                                NCores = 2L,
                                Name = paste0(sample(c(0:9, letters), 10), collapse = ""),
                                Verbose = TRUE) {

  if (! is.list(CountsList)) CountsList <- list(CountsList)

  # filter counts
  filter_mat <- sapply(CountsList, function(x) rowMeans(x >= MinCount) >= PropMinCount)
  filter <- apply(filter_mat, 1, any)
  ercc_filter <- apply(filter_mat[grepl("^ERCC", names(filter)), ], 1, all)
  filter[grepl("^ERCC", names(filter))] <- ercc_filter

  print(c(n_trascripts = sum(filter),
          proportion_transcripts_kept = sum(filter)/length(filter),
          n_ercc = sum(ercc_filter),
          proportion_ercc_kept = sum(ercc_filter)/length(ercc_filter)))

  # remove spike ins that did not pass filter
  SpikeInfo <- SpikeInfo[SpikeID %in% grep("^ERCC", names(filter)[filter], value = TRUE), ]
  SpikeInfo <- as.data.frame(SpikeInfo)

  # format for BASiCS
  dat_list <- lapply(CountsList, function(x)
    suppressMessages(BASiCS::newBASiCS_Data(Counts = x[filter, ],
                                            Tech = grepl("^ERCC", rownames(x[filter, ])),
                                            SpikeInfo = SpikeInfo)))

  if (Verbose) {
    if (! dir.exists("./logs")) dir.create("./logs/")
  }

  # run MCMC
  raw_chain_list <- parallel::mclapply(1:length(dat_list), function(i) {

    if (Verbose) {
      temps <- format(Sys.time(), "%Y_%m_%d_%Hh%M")
      output_filename <- file(paste0("./logs/basics_", Name, temps, "_chain", i, ".out"), open = "wt")
      sink(file = output_filename, type = "output")
    }

    chain <- BASiCS::BASiCS_MCMC(Data = dat_list[[i]],
                                 N = N,
                                 Thin = Thin,
                                 Burn = Burn,
                                 Regression = Regression,
                                 WithSpikes = TRUE,
                                 PrintProgress = TRUE)

    if (Verbose) {
      sink()
      close(output_filename)
    }

    return(chain)
  }, mc.cores = NCores)

  # format and return results
  names(raw_chain_list) <- names(dat_list)
  return(raw_chain_list)
}
