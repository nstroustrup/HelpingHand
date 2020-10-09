#' Cut a Tree into Groups of Data
#'
#' Does the same as \link[stats]{cutree}, but the numbering
#' matches the order of clusters when plotting the tree
#' produced by \link[stats]{hclust}
#'
#' @inheritParams stats::cutree
#'
#' @export
#'
cutree2 <- function(tree, k = NULL, h = NULL) {

  reorderCluster <- function(x, tree) {
    o <- x[match(tree$labels[tree$order], names(x))]
    o <- unique(o)
    y <- numeric(length(x))
    names(y) <- names(x)
    for (i in 1:max(x)) y[x == o[i]] <- i
    y
  }

  groups <- stats::cutree(tree = tree, k = k, h = h)

  if (! is.matrix(groups)) {
    groups2 <- reorderCluster(groups, tree)
  } else {
    groups2 <- apply(groups, 2, reorderCluster, tree = tree)
  }
  groups2

}
