#' @export
normalizeUsingSpikeins <- function(dat, pseudo_count = 0.5) {

  # normalize without offset
  results <- lapply(dat, function(sce) {


    # subset of counts
    bio_x <- counts(sce)
    tech_x <- assay(altExp(sce))

    x <- rbind(bio_x, tech_x)

    tech <- c(rep_len(FALSE, nrow(bio_x)), rep_len(TRUE, nrow(tech_x)))
    bio <- ! tech
    input <- metadata(sce)$SpikeInput[, 2L]

    rm(bio_x, tech_x)

    # technical factor - nu
    counts_tech <- as.data.table(cbind(input = input, x[tech, ]))
    counts_tech <- melt(counts_tech, id.vars = "input", variable.name = "sample", value.name = "count")
    nu <- sapply(colnames(x), function(s) {
      fit <- MASS::glm.nb(formula = count ~ log(input), data = counts_tech[sample == s])
      exp(coef(fit)["(Intercept)"])
    })
    names(nu) <- colnames(x)
    # nu <- counts_tech[count > spikein_mincount, exp(median(log(count) - log(input))), by = sample]$V1

    # size factor - phi
    logprenorm <- log(sapply(1:ncol(x), function(i) (x[bio, i]+pseudo_count) / nu[i]))
    ref <- rowMedians(logprenorm)
    phi <- apply(logprenorm, 2, function(x) exp(median(x - ref)))
    phi <- length(phi) * phi/sum(phi)
    names(phi) <- colnames(x)

    # median mean gene expression
    medexp <- median(rowMedians(sapply(1:ncol(x), function(i) (x[bio, i]) / (phi[i] * nu[i]))))

    list(nu = nu, phi = phi, medexp = medexp)
  })
  names(results) <- names(dat)

  # compute offset
  refmedexp <- results[[1]]$medexp
  offset <- sapply(results, function(x) x$medexp) / refmedexp

  # normalize with offset
  results <- lapply(1:length(dat), function(i) {

    sce <- dat[[i]]

    # subset of counts
    bio_x <- counts(sce)
    tech_x <- assay(altExp(sce))
    x <- rbind(bio_x, tech_x)

    tech <- c(rep_len(FALSE, nrow(bio_x)), rep_len(TRUE, nrow(tech_x)))
    bio <- ! tech

    rm(bio_x, tech_x)

    o <- offset[i]
    phi <- results[[i]]$phi * o
    nu <- results[[i]]$nu

    # normalized counts with offset
    y <- sapply(1:ncol(x), function(i) {
      c(x[bio, i] / (phi[i] * nu[i]), x[tech, i] / nu[i])
    })
    y <- y[match(rownames(x), rownames(y)), ]
    rownames(y) <- rownames(x)
    colnames(y) <- colnames(x)
    y
    list(y = y, nu = nu, phi = phi, offset = o)
  })
  names(results) <- names(dat)

  # format
  y <- Reduce(cbind, lapply(results, `[[`, "y"))
  nu <- Reduce(c, lapply(results, `[[`, "nu"))
  phi <- Reduce(c, lapply(results, `[[`, "phi"))
  o <- Reduce(c, lapply(results, `[[`, "offset"))

  list(y = y, nu = nu, phi = phi, offset = o)

}
