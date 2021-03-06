% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BASiCS_RegressionDE.R
\name{BASiCS_RegressionDE}
\alias{BASiCS_RegressionDE}
\title{Detection of genes with changes in expression using linear regression}
\usage{
BASiCS_RegressionDE(
  Chains,
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
  ...
)
}
\arguments{
\item{Chains}{named list of object of class \code{\linkS4class{BASiCS_Chain}}.
No offset is implemented here. Please use \code{BASiCS_CorrectOffset} before.
.}

\item{EFDR_M}{Target for expected false discovery rate related to
the comparison of means. If \code{EFDR_M = NULL}, EFDR calibration is not
performed and the posterior probability threshold is set equal to 
\code{ProbThresholdM}. Default \code{EFDR_M = 0.05}.}

\item{EFDR_D}{Target for expected false discovery rate related to
the comparison of dispersions. If \code{EFDR_D = NULL}, EFDR calibration is 
not performed and the posterior probability threshold is set equal to 
\code{ProbThresholdD}.Default \code{EFDR_D = 0.05}.}

\item{EFDR_R}{Target for expected false discovery rate related to
the comparison of residual over-dispersions. If \code{EFDR_R = NULL}, EFDR 
calibration is not performed and the posterior probability threshold is set 
equal to \code{ProbThresholdR}.Default \code{EFDR_D = 0.05}.}

\item{EpsilonM}{Minimum fold change tolerance threshold for detecting
changes in overall expression (must be a positive real number).
Default value: \code{EpsilonM = log2(1.5)} (i.e. 50\% increase).}

\item{EpsilonD}{Minimum fold change tolerance threshold for detecting
changes in biological over-dispersion (must be a positive real number).
Default value: \code{EpsilonM = log2(1.5)} (i.e. 50\% increase).}

\item{EpsilonR}{Minimum distance threshold for detecting
changes in residual over-dispersion (must be a positive real number).
Default value: \code{EpsilonR= log2(1.5)/log2(exp(1))} (i.e. 50\% increase).}

\item{ProbThresholdM}{Optional parameter. Probability threshold for detecting
changes in overall expression (must be a positive value, between 0 and 1). 
If \code{EFDR_M = NULL}, the posterior probability threshold for the 
differential mean expression test will be set to \code{ProbThresholdM}. If
a value for \code{EFDR_M} is provided, the posterior probability threshold
is chosen to achieve an EFDR equal to \code{EFDR_M} and \code{ProbThresholdM} 
defines a minimum probability threshold for this calibration (this avoids low 
values of \code{ProbThresholdM} to be chosen by the EFDR calibration. 
Default value \code{ProbThresholdM = 2/3}, i.e. the probability of observing 
a log2-FC above \code{EpsilonM} must be at least twice the probality of 
observing the complementary event (log2-FC below \code{EpsilonM}).}

\item{ProbThresholdD}{Optional parameter. Probability threshold for detecting
changes in cell-to-cell biological over-dispersion (must be a positive value, 
between 0 and 1). Same usage as \code{ProbThresholdM}, depending on the value
provided for \code{EFDR_D}. Default value \code{ProbThresholdD = 2/3}.}

\item{ProbThresholdR}{Optional parameter. Probability threshold for detecting
changes in residual over-dispersion (must be a positive value, between 0 and 
1). Same usage as \code{ProbThresholdM}, depending on the value provided for 
\code{EFDR_R}. Default value \code{ProbThresholdR = 2/3}.}

\item{OrderVariable}{Ordering variable for output.
Possible values: \code{'GeneIndex'} (default), \code{'GeneName'} and 
\code{'Mu'} (mean expression).}

\item{GenesSelect}{Optional argument to provide a user-defined list
of genes to be considered for the comparison.
Default: \code{GenesSelect = NULL}. When used, this argument must be a vector
of \code{TRUE} (include gene) / \code{FALSE} (exclude gene) indicator,
with the same length as the number of intrinsic genes and following the same
order as how genes are displayed in the table of counts.
This argument is necessary in order to have a meaningful EFDR calibration
when the user decides to exclude some genes from the comparison.}

\item{Parameters}{specifies which parameters should be tested.}

\item{...}{Optional parameters.}

\item{Design}{an object of class \code{data.frame} specifiying design of study.
At least one column must be named \code{Chain} and correspond to \code{names(Chains)}.}

\item{Formula}{an object of class \code{formula} specifying a description of the model to be fitted.}

\item{MultiClass}{is a boolean specifying if all coefficients should be tested}
}
\value{
\code{BASiCS_RegressionDE} returns a list similar to BASiCS_TestDE.
}
\description{
Function to assess changes in expression between two or more groups
of cells (mean and over-dispersion).
}
\seealso{
\link[BASiCS]{BASiCS_TestDE}
}
