#' @name reverselog_trans
#'
#' @title Inverse log scale
#'
#' @description Transform data in the inverse log scale for ggplot2
#'
#' @param base numeric; logarithm base.
#'
#' @seealso \link[ggplot2]{scale_continuous}
#'
#' @examples
#' \donttest{
#' ggplot(trees, aes(Girth, Height)) +
#'   geom_point() +
#'   scale_y_continuous(trans=reverselog_trans())
#'   }
#'
#' @import ggplot2
#' @export
#'
reverselog_trans <- function(base = 10) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                    scales::log_breaks(base = base),
                    domain = c(-Inf, Inf))
}
