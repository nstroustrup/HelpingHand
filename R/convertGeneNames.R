#' @name convertGene
#'
#' @title Convert gene names
#'
#' @description Convert wormbase gene identification to official gene symbol. Uses \code{org.Ce.eg.db} for mapping.
#'
#' @inheritParams AnnotationDbi::mapIds
#'
#' @return character; Converted gene identification to the type specified by \code{column}.
#'
#' @seealso \link[AnnotationDbi]{mapIds}
#'
#' @export
#'
convertGene <- function(keys, column = "SYMBOL", keytype = "WORMBASE", x = org.Ce.eg.db::org.Ce.eg.db) {
  suppressMessages(AnnotationDbi::mapIds(x = org.Ce.eg.db::org.Ce.eg.db,
                                         keys = keys,
                                         column = column,
                                         keytype = keytype))
}
