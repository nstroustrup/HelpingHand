#' @name runGSEA
#'
#' @title Run gene set enrichment analysis for C. elegans
#'
#' @description Wrapper function around \link[clusterProfiler]{gseKEGG} and \link[clusterProfiler]{gseGO}.
#'
#' @inheritParams clusterProfiler::gseKEGG
#' @inheritParams clusterProfiler::gseGO
#' @inheritParams AnnotationDbi::mapIds
#'
#' @return A list of \code{gseaResult} for KEGG and GO annotations.
#'
#' @export
#'
runGSEA <- function(geneList,
                    keytype = "WORMBASE",
                    nPerm = 50000L,
                    pvalueCutoff = 1,
                    use_internal_data = TRUE) {


  # Format geneList
  names(geneList) <- convertGene(keys = names(geneList),
                                 column = "ENTREZID",
                                 keytype = keytype)
  geneList <- geneList[! is.na(geneList)]

  geneList <- tapply(geneList, factor(names(geneList)), function(x) x[which.max(abs(x))])
  nms <- names(geneList)
  geneList <- as.numeric(geneList)
  names(geneList) <- nms

  geneList <- sort(geneList, decreasing = TRUE)

  # Run GSEA
  kegg <- .tryCatchNull(clusterProfiler::gseKEGG(geneList          = geneList,
                                                 use_internal_data = use_internal_data,
                                                 organism          = 'cel',
                                                 keyType           = "ncbi-geneid",
                                                 pvalueCutoff      = pvalueCutoff,
                                                 pAdjustMethod     = "BH",
                                                 nPerm             = nPerm,
                                                 minGSSize         = 10,
                                                 maxGSSize         = 150,
                                                 verbose           = FALSE))

  go_bp <- .tryCatchNull(clusterProfiler::gseGO(geneList      = geneList,
                                                OrgDb         = org.Ce.eg.db::org.Ce.eg.db,
                                                ont           = "BP",
                                                nPerm         = nPerm,
                                                pAdjustMethod = "BH",
                                                minGSSize     = 10,
                                                maxGSSize     = 200,
                                                pvalueCutoff  = pvalueCutoff,
                                                verbose       = FALSE))
  if (! is.null(go_bp)) {
    aux <- .convertEntrez(go_bp@result$core_enrichment)
    go_bp@result[, "core_enrichment_symbol"] <- aux$symbol
    go_bp@result[, "core_enrichment_wormbase"] <- aux$wormbase
    go_bp@result[, "Description"] <- .formatNames(go_bp@result[, "Description"])
  }

  go_mf <- .tryCatchNull(clusterProfiler::gseGO(geneList      = geneList,
                                                OrgDb         = org.Ce.eg.db::org.Ce.eg.db,
                                                ont           = "MF",
                                                nPerm         = nPerm,
                                                pAdjustMethod = "BH",
                                                minGSSize     = 10,
                                                maxGSSize     = 200,
                                                pvalueCutoff  = pvalueCutoff,
                                                verbose       = FALSE))

  if (! is.null(go_mf)) {
    aux <- .convertEntrez(go_mf@result$core_enrichment)
    go_mf@result[, "core_enrichment_symbol"] <- aux$symbol
    go_mf@result[, "core_enrichment_wormbase"] <- aux$wormbase
    go_mf@result[, "Description"] <- .formatNames(go_mf@result[, "Description"])
  }

  go_cc <- .tryCatchNull(clusterProfiler::gseGO(geneList      = geneList,
                                                OrgDb         = org.Ce.eg.db::org.Ce.eg.db,
                                                ont           = "CC",
                                                nPerm         = nPerm,
                                                pAdjustMethod = "BH",
                                                minGSSize     = 10,
                                                maxGSSize     = 200,
                                                pvalueCutoff  = pvalueCutoff,
                                                verbose       = FALSE))
  if (! is.null(go_cc)) {

    aux <- .convertEntrez(go_cc@result$core_enrichment)
    go_cc@result[, "core_enrichment_symbol"] <- aux$symbol
    go_cc@result[, "core_enrichment_wormbase"] <- aux$wormbase
    go_cc@result[, "Description"] <- .formatNames(go_cc@result[, "Description"])
  }

  tables <- list(kegg = kegg, go_bp = go_bp, go_cc = go_cc, go_mf = go_mf, nPerm = nPerm)
  tables
}
#' @name plotGSEA
#'
#' @title Plot GSEA results
#'
#' @description Plot results produced by \link{runGSEA}.
#'
#' @param tables List of \code{enrichResult} produced by \link{runGSEA}.
#' @inheritParams clusterProfiler::ridgeplot
#'
#' @return An object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
#'
plotGSEA <- function(tables,
                     showCategory = 50L) {

  kegg <- tables$kegg
  go_bp <- tables$go_bp
  go_cc <- tables$go_cc
  go_mf <- tables$go_mf
  nPerm <- tables$nPerm

  if (! is.null(go_bp)) {
    go_bp_plot <- suppressWarnings(clusterProfiler::ridgeplot(go_bp, showCategory = showCategory) +
      ggplot2::scale_fill_gradient2(trans = "log10", limit = c(1/nPerm, 1), midpoint = log10(0.01))) +
      ggplot2::ggtitle(paste("GO BP")) +
      ggplot2::xlab("Enrichment score") +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(title = element_text(size = 15))
  } else {
    go_bp_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  }

  if (! is.null(go_mf)) {
    go_mf_plot <- suppressWarnings(clusterProfiler::ridgeplot(go_mf, showCategory = showCategory) +
      ggplot2::scale_fill_gradient2(trans = "log10", limit = c(1/nPerm, 1), midpoint = log10(0.01))) +
      ggplot2::ggtitle(paste("GO MF")) +
      ggplot2::xlab("Enrichment score") +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(title = element_text(size = 15))
  } else {
    go_mf_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  }

  if (! is.null(go_cc)) {
    go_cc_plot <- suppressWarnings(clusterProfiler::ridgeplot(go_cc, showCategory = showCategory) +
      ggplot2::scale_fill_gradient2(trans = "log10", limit = c(1/nPerm, 1), midpoint = log10(0.01))) +
      ggplot2::ggtitle(paste("GO CC")) +
      ggplot2::xlab("Enrichment score") +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(title = element_text(size = 15))
  } else {
    go_cc_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  }

  if (! is.null(kegg)) {
    aux <- .convertEntrez(kegg@result$core_enrichment)
    kegg@result[, "core_enrichment_symbol"] <- aux$symbol
    kegg@result[, "core_enrichment_wormbase"] <- aux$wormbase

    kegg_plot <- suppressWarnings(clusterProfiler::ridgeplot(kegg, showCategory = showCategory) +
      ggplot2::scale_fill_gradient2(trans = "log10", limit = c(1/nPerm, 1), midpoint = log10(0.01))) +
      ggplot2::ggtitle(paste("KEGG")) +
      ggplot2::xlab("Enrichment score") +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(title = element_text(size = 15))
  } else {
    kegg_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  }

  plots <- suppressMessages(ggpubr::ggarrange(go_bp_plot, go_mf_plot, go_cc_plot, kegg_plot, nrow = 2L, ncol = 2L))
  plots <- ggpubr::annotate_figure(plots, top = ggpubr::text_grob("Gene set enrichment analysis\n", size = 18))

  plots
}

#' @name runEnrichment
#'
#' @title Run enrichment analysis of a gene set for C. elegans
#'
#' @description Wrapper function around \link[clusterProfiler]{enrichKEGG} and \link[clusterProfiler]{enrichGO}.
#'
#' @inheritParams clusterProfiler::enrichKEGG
#' @inheritParams clusterProfiler::enrichGO
#' @inheritParams AnnotationDbi::mapIds
#'
#' @return A list of \code{enrichResult} for KEGG and GO annotations.
#'
#' @export
#'
runEnrichment <- function(gene, universe = NULL, keytype = "WORMBASE", pvalueCutoff = 1, use_internal_data = TRUE) {

  # Format gene set and universe
  gene <- convertGene(keys = gene,
                      column = "ENTREZID",
                      keytype = keytype)
  gene <- unique(gene[! is.na(gene)])

  universe <- convertGene(keys = universe,
                          column = "ENTREZID",
                          keytype = keytype)
  universe <- unique(universe[! is.na(universe)])

  # Run GSEA
  kegg <- clusterProfiler::enrichKEGG(gene = gene,
                                      use_internal_data = use_internal_data,
                                      organism  = "cel",
                                      universe = universe,
                                      pvalueCutoff = pvalueCutoff,
                                      keyType = "ncbi-geneid",
                                      pAdjustMethod = "BH",
                                      minGSSize = 5,
                                      maxGSSize = 150)
  if (! is.null(kegg)) {
    aux <- .convertEntrez(kegg@result$geneID)
    kegg@result[, "geneID_symbol"] <- aux$symbol
    kegg@result[, "geneID_wormbase"] <- aux$wormbase
  }

  go_bp <- clusterProfiler::enrichGO(gene          = gene,
                                     OrgDb         = org.Ce.eg.db::org.Ce.eg.db,
                                     ont           = "BP",
                                     universe      = universe,
                                     pvalueCutoff  = pvalueCutoff,
                                     pAdjustMethod = "BH",
                                     minGSSize     = 10,
                                     maxGSSize     = 200)
  if (! is.null(go_bp)) {
    aux <- .convertEntrez(go_bp@result$geneID)
    go_bp@result[, "geneID_symbol"] <- aux$symbol
    go_bp@result[, "geneID_wormbase"] <- aux$wormbase
    go_bp@result[, "Description"] <- .formatNames(go_bp@result[, "Description"])
  }


  go_mf <- clusterProfiler::enrichGO(gene          = gene,
                                     OrgDb         = org.Ce.eg.db::org.Ce.eg.db,
                                     ont           = "MF",
                                     universe      = universe,
                                     pvalueCutoff  = pvalueCutoff,
                                     pAdjustMethod = "BH",
                                     minGSSize     = 10,
                                     maxGSSize     = 200)
  if (! is.null(go_mf)) {
    aux <- .convertEntrez(go_mf@result$geneID)
    go_mf@result[, "geneID_symbol"] <- aux$symbol
    go_mf@result[, "geneID_wormbase"] <- aux$wormbase
    go_mf@result[, "Description"] <- .formatNames(go_mf@result[, "Description"])
  }

  go_cc <- clusterProfiler::enrichGO(gene          = gene,
                                     OrgDb         = org.Ce.eg.db::org.Ce.eg.db,
                                     ont           = "CC",
                                     universe      = universe,
                                     pvalueCutoff  = pvalueCutoff,
                                     pAdjustMethod = "BH",
                                     minGSSize     = 10,
                                     maxGSSize     = 200)

  if (! is.null(go_cc)) {
    aux <- .convertEntrez(go_cc@result$geneID)
    go_cc@result[, "geneID_symbol"] <- aux$symbol
    go_cc@result[, "geneID_wormbase"] <- aux$wormbase
    go_cc@result[, "Description"] <- .formatNames(go_cc@result[, "Description"])
  }

  tables <- list(kegg = kegg, go_bp = go_bp, go_cc = go_cc, go_mf = go_mf)
  tables
}

#' @name plotEnrichment
#'
#' @title Plot enrichment results
#'
#' @description Plot results produced by \link{runEnrichment}.
#'
#' @param tables List of \code{enrichResult} produced by \link{runEnrichment}.
#' @inheritParams clusterProfiler::dotplot
#'
#' @return An object of class ggarrange, which is a ggplot or a list of ggplot.
#'
#' @export
#'
plotEnrichment <- function(tables,
                           showCategory = 50L) {

  kegg <- tables$kegg
  go_bp <- tables$go_bp
  go_cc <- tables$go_cc
  go_mf <- tables$go_mf

  if (! is.null(kegg)) {
    kegg_plot <- clusterProfiler::dotplot(kegg, showCategory = showCategory, orderBy = "GeneRatio") +
      ggplot2::ggtitle(paste("KEGG")) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(title = element_text(size = 15))
  } else {
    kegg_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  }

  if (! is.null(go_bp)) {
    go_bp_plot <- clusterProfiler::dotplot(go_bp, showCategory = showCategory, orderBy = "GeneRatio") +
      ggplot2::ggtitle(paste("GO BP")) +
      ggplot2::xlab("Enrichment score") +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(title = element_text(size = 15))
  } else {
    go_bp_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  }

  if (! is.null(go_mf)) {
    go_mf_plot <- clusterProfiler::dotplot(go_mf, showCategory = showCategory, orderBy = "GeneRatio") +
      ggplot2::ggtitle(paste("GO MF")) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(title = element_text(size = 15))
  } else {
    go_mf_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  }
  if (! is.null(go_cc)) {
    go_cc_plot <- clusterProfiler::dotplot(go_cc, showCategory = showCategory, orderBy = "GeneRatio") +
      ggplot2::ggtitle(paste("GO CC")) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::theme(title = element_text(size = 15))
  } else {
    go_cc_plot <- ggplot2::ggplot() + ggplot2::theme_void()
  }

  plots <- suppressMessages(ggpubr::ggarrange(go_bp_plot,
                                              go_mf_plot,
                                              go_cc_plot,
                                              kegg_plot,
                                              nrow = 2L,
                                              ncol = 2L))
  plots
}

.formatNames <- function(x) {
  x %>%
    stringr::str_sub(start = 0L, 65L) %>%
    make.names(unique = TRUE) %>%
    gsub(pattern = "[.]", replacement = " ") %>%
    gsub(pattern = "[[:space:]]{2,}", replacement = "")
}


.convertEntrez <- function(x) {
  if (length(x) == 0 | is.null(x)) return(x)
  s <- stringr::str_split(x, "/")
  map <- data.table::data.table(ENTREZ = sort(unique(unlist(s))))
  map[, SYMBOL := convertGene(keys = ENTREZ, keytype = "ENTREZID",  column = "SYMBOL")]
  map[, WORMBASE := convertGene(keys = ENTREZ, keytype = "ENTREZID",  column= "WORMBASE")]
  setkey(map, ENTREZ)
  list(symbol = sapply(s, function(e) paste(map[e]$SYMBOL, collapse = "/")),
       wormbase =  sapply(s, function(e) paste(map[e]$WORMBASE, collapse = "/")))
}

.tryCatchNull <- function(x) {
  tryCatch(x, error = function(e) return(NULL))
}
