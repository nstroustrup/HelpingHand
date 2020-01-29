#' @name BASiCS_RegressionDE
#'
#' @title Detection of genes with changes in expression using linear regression
#'
#' @description Function to assess changes in expression between two or more groups
#' of cells (mean and over-dispersion).
#'
#' @inheritParams BASiCS::BASiCS_TestDE
#'
#' @param Chains named list of object of class \code{\linkS4class{BASiCS_Chain}}.
#' No offset is implemented here. Please use \code{BASiCS_CorrectOffset} before.
#' .
#' @param Design an object of class \code{data.frame} specifiying design of study.
#' At least one column must be named \code{Chain} and correspond to \code{names(Chains)}.
#'
#' @param Formula an object of class \code{formula} specifying a description of the model to be fitted.
#'
#' @param MultiClass is a boolean specifying if all coefficients should be tested
#'
#' @param Parameters specifies which parameters should be tested.
#'
#' @return \code{BASiCS_RegressionDE} returns a list similar to BASiCS_TestDE.
#'
#' @seealso \link[BASiCS]{BASiCS_TestDE}
#'
#' @export
#'
BASiCS_RegressionDE <- function(Chains,
                                # Design,
                                # Formula,
                                ModelMatrix,
                                EFDR_M = 0.05,
                                EFDR_D = 0.05,
                                EFDR_R = 0.05,
                                EpsilonM = log2(1.5),
                                EpsilonD = log2(1.5),
                                EpsilonR = log2(1.5)/log2(exp(1)),
                                ProbThresholdM = 2/3,
                                ProbThresholdD = 2/3,
                                ProbThresholdR = 2/3,
                                OrderVariable = "GeneIndex",
                                GenesSelect = NULL,
                                Classes = as.list(colnames(ModelMatrix)),
                                Parameters = c("mu", "delta", "epsilon"),
                                ...) {


  nChains <- length(Chains)
  GeneName <- colnames(Chains[[1]]@parameters$mu)
  nGenes <- length(GeneName)
  GeneIndex <- seq_len(nGenes)
  nIters <- nrow(Chains[[1]]@parameters$mu)
  nSamples <- sapply(Chains, function(x) ncol(x@parameters$phi))
  if (any(! Parameters %in% c("mu", "delta", "epsilon"))) stop("Parameters must be 'mu'', 'delta'' or 'epsilon'.")

  ResDispExists <- all(sapply(Chains, function(x) !is.null(x@parameters$epsilon))) # mean-variance regression was performed
  RegMean <- "mu" %in% Parameters
  RegDisp <- "delta" %in% Parameters
  RegResDisp <- "epsilon" %in% Parameters & ResDispExists

  # all GeneNames are identical
  if (! all(sapply(Chains[2:nChains], function(x) identical(GeneName, colnames(x@parameters$mu))))) {
    stop("Error: GeneNames do not match between chains.")
  }

  # all chains have the same number of iterations
  if (! all(sapply(Chains[2:nChains], function(x) identical(nIters, nrow(x@parameters$mu))))) {
    stop("Error: Chains do not have the name number of iterations.")
  }

  if (any(! unlist(Classes) %in% colnames(ModelMatrix))) stop("Invalid 'Classes': not in colnames(ModelMatrix).")
  # transform design into model matrix
  # if (missing(ModelMatrix)) {
  # Design <- Design[match(Chain, names(Chains))]
  # ModelMatrix <- model.matrix(Formula, Design)
  # intercept <- colnames(ModelMatrix) == "(Intercept)"
  # if (sum(! intercept) == 0) stop("Can not work on a model only with an intercept.")
  # colnames(ModelMatrix)[intercept] <- "Intercept"
  # colnames(ModelMatrix) <- make.names(colnames(ModelMatrix))
  # }

  # place parameters into array
  if (RegMean) AllMean <- simplify2array(lapply(Chains, function(x) x@parameters$mu))
  if (RegDisp) AllDisp <- simplify2array(lapply(Chains, function(x) x@parameters$delta))
  if (RegResDisp) AllResDisp <- simplify2array(lapply(Chains, function(x) x@parameters$epsilon))

  # compute median parameters for each chain
  if (RegMean) {
    MeanOverall <- vapply(Chains, function(x) matrixStats::colMedians(x@parameters$mu), FUN.VALUE = numeric(nGenes))
    rownames(MeanOverall) <- GeneName
  }
  if (RegDisp) {
    DispOverall <- vapply(Chains, function(x) matrixStats::colMedians(x@parameters$delta), FUN.VALUE = numeric(nGenes))
    rownames(DispOverall) <- GeneName
  }
  if (RegResDisp) {
    ResDispOverall <- vapply(Chains, function(x) matrixStats::colMedians(x@parameters$epsilon), FUN.VALUE = numeric(nGenes))
    rownames(ResDispOverall) <- GeneName
  }

  # weight median according to sample size
  if (RegMean) MeanOverall <- apply(MeanOverall, 1, function(x) sum(nSamples * x)) / sum(nSamples)
  if (RegDisp) DispOverall <- apply(DispOverall, 1, function(x) sum(nSamples * x)) / sum(nSamples)
  if (RegResDisp) ResDispOverall <- apply(ResDispOverall, 1, function(x) sum(nSamples * x)) / sum(nSamples)

  # convert to log2
  if (RegMean) AllMean <- log2(AllMean)
  if (RegDisp) AllDisp <- log2(AllDisp)

  # perform linear regression on mu
  if (RegMean) {
    message("Performing linear regressions on mean\n")
    MeanCoefficients <- simplify2array(pbapply::pblapply(GeneIndex, function(gene) {
      sapply(1L:nIters, function(iteration) {
        .lm.fit(ModelMatrix, AllMean[iteration, gene, ])$coefficients
      })
    }))
    dimnames(MeanCoefficients) <- list(colnames(ModelMatrix), NULL, GeneName)

    MedianMeanCoefficient <- apply(MeanCoefficients, 1, function(x) matrixStats::colMedians(x))
    rownames(MedianMeanCoefficient) <- GeneName
    colnames(MedianMeanCoefficient) <- colnames(ModelMatrix)

    SdMeanCoefficient <- apply(MeanCoefficients, 1, function(x) matrixStats::colSds(x))
    rownames(SdMeanCoefficient) <- GeneName
    colnames(SdMeanCoefficient) <- colnames(ModelMatrix)

  } else {
    MeanCoefficients <- NULL
  }

  # perform linear regression on delta
  if (RegDisp) {
    message("Performing linear regressions on dispersion\n")
    DispCoefficients <- simplify2array(pbapply::pblapply(GeneIndex, function(gene) {
      sapply(1L:nIters, function(iteration) {
        .lm.fit(ModelMatrix, AllDisp[iteration, gene, ])$coefficients
      })
    }))
    dimnames(DispCoefficients) <- list(colnames(ModelMatrix), NULL, GeneName)

    MedianDispCoefficient <- apply(DispCoefficients, 1, function(x) matrixStats::colMedians(x))
    rownames(MedianDispCoefficient) <- GeneName
    colnames(MedianDispCoefficient) <- colnames(ModelMatrix)

    SdDispCoefficient <- apply(DispCoefficients, 1, function(x) matrixStats::colSds(x))
    rownames(SdDispCoefficient) <- GeneName
    colnames(SdDispCoefficient) <- colnames(ModelMatrix)

  } else {
    DispCoefficients <- NULL
  }

  # perform linear regression on epsilon
  if (RegResDisp) {
    message("Performing linear regressions on residual dispersion\n")
    ResDispCoefficients <- simplify2array(pbapply::pblapply(GeneIndex, function(gene) {
      sapply(1L:nIters, function(iteration) {
        .lm.fit(ModelMatrix, AllResDisp[iteration, gene, ])$coefficients
      })
    }))
    dimnames(ResDispCoefficients) <- list(colnames(ModelMatrix), NULL, GeneName)

    MedianResDispCoefficient <- apply(ResDispCoefficients, 1, function(x) matrixStats::colMedians(x))
    rownames(MedianResDispCoefficient) <- GeneName
    colnames(MedianResDispCoefficient) <- colnames(ModelMatrix)

    SdResDispCoefficient <- apply(ResDispCoefficients, 1, function(x) matrixStats::colSds(x))
    rownames(SdResDispCoefficient) <- GeneName
    colnames(SdResDispCoefficient) <- colnames(ModelMatrix)

  } else {
    ResDispCoefficients <- NULL
  }

  # compute probabilities and find thresholds
  ClassNames <- sapply(Classes, paste0, collapse="_")
  if (RegMean) {
    AuxMean <- lapply(Classes, function(CoeffNames) {

      HiddenThresholdSearchRegressionDE(Chain = MeanCoefficients[CoeffNames, , ],
                                        Epsilon = EpsilonM,
                                        ProbThreshold = ProbThresholdM,
                                        GenesSelect = GenesSelect,
                                        EFDR = EFDR_M,
                                        Task = paste0("Differential mean (", paste(CoeffNames, collapse = " "), ")"),
                                        Suffix = "M")
    })
    names(AuxMean) <- ClassNames
  }

  if (RegDisp) {
    AuxDisp <- lapply(Classes, function(CoeffNames) {

      HiddenThresholdSearchRegressionDE(Chain = DispCoefficients[CoeffNames, , ],
                                        Epsilon = EpsilonD,
                                        ProbThreshold = ProbThresholdD,
                                        GenesSelect = GenesSelect,
                                        EFDR = EFDR_D,
                                        Task = paste0("Differential dispersion (", paste(CoeffNames, collapse = " "), ")"),
                                        Suffix = "D")
    })
    names(AuxDisp) <- ClassNames
  }

  if (RegResDisp) {
    AuxResDisp <- lapply(Classes, function(CoeffNames) {

      HiddenThresholdSearchRegressionDE(Chain = ResDispCoefficients[CoeffNames, , ],
                                        Epsilon = EpsilonR,
                                        ProbThreshold = ProbThresholdR,
                                        GenesSelect = GenesSelect,
                                        EFDR = EFDR_R,
                                        Task = paste0("Differential residual dispersion (", paste(CoeffNames, collapse = " "), ")"),
                                        Suffix = "R")
    })
    names(AuxResDisp) <- ClassNames
  }

  # summary of thresholds
  if (RegMean) {
    DiffMeanSummary <- sapply(AuxMean, function(x) x$OptThreshold)
    rownames(DiffMeanSummary) <- c("ProbThreshold", "EFDR", "EFNR")
  } else {
    DiffMeanSummary <- NULL
  }

  if (RegDisp) {
    DiffDispSummary <- sapply(AuxDisp, function(x) x$OptThreshold)
    rownames(DiffDispSummary) <- c("ProbThreshold", "EFDR", "EFNR")
  } else {
    DiffDispSummary <- NULL
  }

  if (RegResDisp) {
    DiffResDispSummary <- sapply(AuxResDisp, function(x) x$OptThreshold)
    rownames(DiffResDispSummary) <- c("ProbThreshold", "EFDR", "EFNR")
  } else {
    DiffResDispSummary <- NULL
  }

  # tables
  if (RegMean) {
    TableMean <- data.table::data.table(GeneName = GeneName, MeanOverall = MeanOverall)
    TableMean <- cbind(TableMean, .HiddenFormatAux(AuxMean,
                                                   MedianMeanCoefficient,
                                                   SdMeanCoefficient,
                                                   GenesSelect))
  } else {
    TableMean <- NULL
  }

  if (RegDisp) {
    TableDisp <- data.table::data.table(GeneName = GeneName, MeanOverall = MeanOverall, DispOverall = DispOverall)
    TableDisp <- cbind(TableDisp, .HiddenFormatAux(AuxDisp,
                                                   MedianDispCoefficient,
                                                   SdDispCoefficient,
                                                   GenesSelect))
  } else {
    TableDisp <- NULL
  }

  if (RegResDisp) {
    TableResDisp <- data.table::data.table(GeneName = GeneName, MeanOverall = MeanOverall, ResDispOverall = ResDispOverall)
    TableResDisp <- cbind(TableResDisp, .HiddenFormatAux(AuxResDisp,
                                                         MedianResDispCoefficient,
                                                         SdResDispCoefficient,
                                                         GenesSelect))
  } else {
    TableResDisp <- NULL
  }


  list(TableMean = TableMean,
       TableDisp = TableDisp,
       TableResDisp = TableResDisp,
       MeanCoefficients = MeanCoefficients,
       DispCoefficients = DispCoefficients,
       ResDispCoefficients = ResDispCoefficients,
       DiffMeanSummary = DiffMeanSummary,
       DiffDispSummary = DiffDispSummary,
       DiffResDispSummary = DiffResDispSummary)

}

.HiddenFormatAux <- function(Aux,
                             MedianCoefficient,
                             SdCoefficient,
                             GenesSelect) {

  ClassNames <- names(Aux)
  Table <- lapply(ClassNames, function(Class) {
    A <- Aux[[Class]]
    Prob <- A$Prob
    OptThreshold <- A$OptThreshold

    DE <- which(Prob > OptThreshold[1])
    ResultDiff <- rep("NoDiff", length(Prob))
    ResultDiff[DE] <- "Diff"
    if (!is.null(GenesSelect)) {
      ResultDiff[!GenesSelect] <- "ExcludedByUser"
    }

    Table <- data.table::data.table(ProbDiffMean = Prob,
                                    ResultDiff = ResultDiff)
    colnames(Table) <- paste(Class, colnames(Table), sep = "_")
    Table
  })

  MedianCoefficientTable <- data.table::as.data.table(MedianCoefficient)
  colnames(MedianCoefficientTable) <- paste0(colnames(MedianCoefficientTable), "_Coefficient")

  SdCoefficientTable <- data.table::as.data.table(SdCoefficient)
  colnames(SdCoefficientTable) <- paste0(colnames(SdCoefficientTable), "_SD")

  Table <- cbind(MedianCoefficientTable, cbind(SdCoefficientTable, Reduce(cbind, Table)))

  return(Table)
}

HiddenTailProbRegressionDE <- function(Chain,
                                       Epsilon) {

  if (Epsilon > 0) {
    Aux <- abs(Chain) > Epsilon
    if (length(dim(Aux)) > 2) Aux <- Reduce(function(x, y) x & y, lapply(1L:dim(Aux)[1L], function(i) Aux[i, , ]))
    Prob <- matrixStats::colMeans2(Aux)

  } else {
    Aux <- Chain > 0
    if (length(dim(Aux)) > 2) Aux <- Reduce(function(x, y) x & y, lapply(1L:dim(Aux)[1L], function(i) Aux[i, , ]))
    Prob_aux <- matrixStats::colMeans2(Aux)
    Prob <- 2 * pmax(Prob_aux, 1 - Prob_aux) - 1
  }

  return(Prob)
}

HiddenThresholdSearchRegressionDE <- function(Chain,
                                              Epsilon,
                                              ProbThreshold,
                                              GenesSelect,
                                              EFDR,
                                              Task,
                                              Suffix) {

  # Calculating posterior probabilities
  Prob <- HiddenTailProbRegressionDE(Chain, Epsilon)

  # Posterior probability threshold search
  if (!is.null(EFDR)) {

    ProbThresholds <- seq(0.5, 0.9995, by = 0.00025)

    if (is.null(GenesSelect)) {
      EFDRgrid <- vapply(ProbThresholds, FUN = HiddenEFDR,
                         FUN.VALUE = 1, Prob = Prob)
      EFNRgrid <- vapply(ProbThresholds, FUN = HiddenEFNR,
                         FUN.VALUE = 1, Prob = Prob)
    }
    else {
      EFDRgrid <- vapply(ProbThresholds, FUN = HiddenEFDR,
                         FUN.VALUE = 1, Prob = Prob[GenesSelect])
      EFNRgrid <- vapply(ProbThresholds, FUN = HiddenEFNR,
                         FUN.VALUE = 1, Prob = Prob[GenesSelect])
    }

    above <- abs(EFDRgrid - EFDR)

    if (sum(!is.na(above)) > 0) {

      # Search EFDR closest to the desired value
      EFDRopt <- EFDRgrid[above == min(above, na.rm = TRUE) & !is.na(above)]
      # If multiple threholds lead to same EFDR, choose the one with lowest EFNR
      EFNRopt <- EFNRgrid[EFDRgrid == mean(EFDRopt) & !is.na(EFDRgrid)]
      if (sum(!is.na(EFNRopt)) > 0) {
        optimal <- which(EFDRgrid == mean(EFDRopt) & EFNRgrid == mean(EFNRopt))
      }
      else {
        optimal <- which(EFDRgrid == mean(EFDRopt))
      }
      # Quick fix for EFDR/EFNR ties; possibly not an issue in real datasets
      optimal <- median(round(median(optimal)))

      # If calibrated threshold is above the minimum required probability
      if(ProbThresholds[optimal] >= ProbThreshold) {

        OptThreshold <- c(ProbThresholds[optimal],
                          EFDRgrid[optimal], EFNRgrid[optimal])
        if (abs(OptThreshold[2] - EFDR) > 0.025) {
          message("For ", Task, " task:\n",
                  "It is not possible to find a probability threshold (>0.5) \n",
                  "that achieves the desired EFDR level (+-0.025). \n",
                  "The output below reflects the closest possible value. \n")
        }
      }
      # If calibrated threshold is below the minimum required probability
      else {
        if (is.null(GenesSelect)) {
          EFDRgrid <- HiddenEFDR(ProbThreshold, Prob)
          EFNRgrid <- HiddenEFNR(ProbThreshold, Prob)
        }
        else {
          EFDRgrid <- HiddenEFDR(ProbThreshold, Prob[GenesSelect])
          EFNRgrid <- HiddenEFNR(ProbThreshold, Prob[GenesSelect])
        }
        OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])
        message("For ", Task, " task:\n",
                "the posterior probability threshold chosen via EFDR calibration is too low.",
                "Probability threshold automatically set equal to 'ProbThreshold", Suffix, "'.")
      }
    }
    else {
      message("EFDR calibration failed for ", Task, " task. \n",
              "Probability threshold automatically set equal to 0.90 \n")
      OptThreshold <- c(0.9, NA, NA)
    }
  }
  else {
    # When a posterior probability threshold has been set a priori
    if (is.null(GenesSelect)) {
      EFDRgrid <- HiddenEFDR(ProbThreshold, Prob)
      EFNRgrid <- HiddenEFNR(ProbThreshold, Prob)
    }
    else {
      EFDRgrid <- HiddenEFDR(ProbThreshold, Prob[GenesSelect])
      EFNRgrid <- HiddenEFNR(ProbThreshold, Prob[GenesSelect])
    }
    OptThreshold <- c(ProbThreshold, EFDRgrid[1], EFNRgrid[1])
  }

  list(Prob = Prob, OptThreshold = OptThreshold,
       EFDRgrid = EFDRgrid, EFNRgrid = EFNRgrid)
}

HiddenEFDR <- function(EviThreshold, Prob) {
  i <- Prob > EviThreshold
  sum((1 - Prob) * i) / sum(i)
  # return(sum((1 - Prob) * I(Prob > EviThreshold)) / sum(I(Prob > EviThreshold)))
}

HiddenEFNR <- function(EviThreshold, Prob) {
  i <- Prob > EviThreshold
  sum(Prob * i) / sum(i)
  # return(sum(Prob * I(EviThreshold >= Prob)) / sum(I(EviThreshold >= Prob)))
}

