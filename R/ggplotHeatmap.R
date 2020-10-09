#' Correlation plot with dendograms
#'
#' Plots a correlation plot with dendograms using ggplot2
#'
#' @param x data matrix
#' @param hc hierarchical clustering (output of \link{hclust})
#' @param k number of clusters
#' @param dendogram_colors cluster colors in dendogram (length should be k)
#'
#' @examples
#' hc <- hclust(as.dist(1-cor(trees)))
#' clusteredCorrplot(trees, hc, 2, c("red", "blue"))
#'
#' @export
#'
clusteredCorrplot <- function(x, hc, k = NULL, dendogram_colors = NULL) {

  # plot dendogram
  dend <- as.dendrogram(hc)
  dend <- dendextend::set(dend, "branches_lwd", 0.4)

  if (! is.null(k)) {
    if (is.null(dendogram_colors)) {
      dend <- dendextend::set(dend, "branches_k_color", k = k)
    } else {
      if (k != length(dendogram_colors))
        stop("Error: length(dendogram_colors) != k")
      dend <- dendextend::set(dend, "branches_k_color", k = k, value = dendogram_colors)
    }
  }

  ggdend <- dendextend::as.ggdend(dend)
  xlims <- c(0.5, nrow(ggdend$labels) + 0.5)

  py <- ggplot(dend, labels = FALSE, na.rm = TRUE) +
    theme_bw() +
    coord_cartesian(expand = FALSE, xlim = xlims) +
    theme(
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.background = element_blank())

  # plot heatmap
  rho <- cor(x)
  rho <- rho[hc$order, ]
  rho <- rho[, hc$order]
  rownames(rho) <- NULL
  colnames(rho) <- NULL

  mlt <- reshape2::melt(rho)
  colnames(mlt) <- c("row", "column", "value")
  mlt <- as.data.table(mlt)

  p <- ggplot2::ggplot(mlt, aes(x = column, y = row, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(name = "Correlation",
                         limits = c(-1, 1),
                         low = "#990000ff", high = "#0b5394ff") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.line = element_line(color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()) +
    coord_cartesian(expand = FALSE)

  plots <- list(p = p,
                px = NULL,
                py = py,
                pr = NULL,
                pc = NULL)

  heatmaply_arrange_plots(
    plots,
    widths = NULL,
    heights = NULL,
    row_dend_left = FALSE)
}

# adapted from heatmaply
heatmaply_arrange_plots <- function(
  plots,
  widths = NULL,
  heights = NULL,
  row_dend_left = FALSE) {


  `%||%` <- function(a, b) {
    if (!is.null(a)) {
      a
    } else {
      b
    }
  }

  ggplot_empty <- function() {
    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::theme(plot.margin = unit(c(0, 0, 0, 0), "npc"))
  }


  default_dims <- function(px, pr) {
    if (!is.null(px)) {
      if (is.null(pr)) {
        widths <- c(0.8, 0.2)
      } else {
        widths <- c(0.7, 0.1, 0.2)
      }
    } else {
      if (is.null(pr)) {
        widths <- 1
      } else {
        widths <- c(0.9, 0.1)
      }
    }
    widths
  }

  plots <- plots[!sapply(plots, is.null)]
  if (!row_dend_left) {
    plots$p <- plots$p + theme(legend.position = "left")
  }
  plots <- lapply(plots, function(x) x +
                    theme(plot.margin = unit(c(0, 0, 0, 0), "npc")))
  plots$py <- plots$py +
    theme(plot.margin = unit(c(1, 1, 0, 0), "cm"))

  plots$p <- plots$p +
    theme(plot.margin = unit(c(0, 1, 1, 1), "cm"))

  column_list <- list(plots$py, plots$pc, plots$p)
  ind_null_col <- sapply(column_list, is.null)

  row1_list <- list(plots$py, ggplot_empty(), ggplot_empty())
  row2_list <- list(plots$pc, ggplot_empty(), ggplot_empty())
  row3_list <- list(plots$p, plots$pr, plots$px)

  if (row_dend_left) {
    row3_list <- rev(row3_list)
    row2_list <- rev(row2_list)
    row1_list <- rev(row1_list)
  }
  plotlist <- c(
    row1_list,
    row2_list,
    row3_list
  )

  nrows <- sum(!ind_null_col)
  ind_remove_col <- rep(ind_null_col, each = length(plotlist) / 3)

  ind_null_row <- sapply(row3_list, is.null)
  ncols <- sum(!ind_null_row)
  ind_remove_row <- rep(ind_null_row, length.out = length(plotlist))
  plotlist <- plotlist[!(ind_remove_row | ind_remove_col)]

  egg::ggarrange(
    plots = plotlist,
    ncol = ncols,
    widths = widths %||% default_dims(plots$px, plots$pr),
    heights = heights %||% rev(default_dims(plots$py, plots$pc))
  )
}
