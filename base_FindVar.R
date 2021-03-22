FastExpMean <- function(x) {
  x <- exp(x) - 1
  rowm <- apply(x, 1, mean)
  rowml <- log(rowm + 1)
  return(rowml)
}
FastLogVMR <- function(x) {
  x <- exp(x) - 1
  rowm <- apply(x, 1, mean)
  nnZero <- apply(x, 1, function(x) sum(x!=0))
  rowv <- apply(x, 1, function(x) (x - mean(x))^2)
  rowv <- (rowv + (ncol(x) - nnZero) * rowm^2) / (ncol(x) - 1)
  rowv <- log(rowv/rowm)
  return(rowv)
}
FastLogVMR <- function(x) {
  x <- exp(x) - 1
  rowm <- apply(x, 1, mean)
  rowv <- apply(x, 1, var)
  rowv <- log(rowv/rowm)
  return(rowv)
}

#' Find variable genes
#'
#' @export
FindVariableFeatures1 <- function(object,
                                 selection.method = c('vst','mean.var.plot','dispersion'),
                                 loess.span = 0.3,
                                 clip.max = 'auto',
                                 num.bin = 20,
                                 binning.method = c('equal_width', 'equal_frequency'),
                                 nfeatures = 2000,
                                 mean.cutoff = c(0.1, 8),
                                 dispersion.cutoff = c(1, Inf),
                                 assay = 'RNA'
) {
  selection.method <- match.arg(selection.method)
  binning.method <- match.arg(binning.method)
  if (selection.method == 'vst') {
    ###1. vst method.
    #mean
    data <- object@assays[[assay]]@counts
    hvf.info <- data.frame(mean = Matrix::rowMeans(x = data))
    #variance
    #hvf.info$variance <- rowVars(as.matrix(data))
    #as.matrix may cause memory problem when used for large data.
    hvf.info$variance <- apply(data, 1, var)
    #variance.expected
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = loess.span
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    #variance.standardized
    sdnor <- function(x, mean.x, var.ex, clip.max){
      nor_x <- (x - mean.x)/var.ex
      #when value in var.ex equal to 0 ,return na. Then convert na to 0.
      nor_x[is.na(nor_x)] <- 0
      nor_x[nor_x > clip.max] <- clip.max
      sdx <- apply(nor_x, 1, var)
      return(sdx)
    }
    if (clip.max == 'auto') {
      clip.max <- sqrt(x = ncol(x = data))
    }
    hvf.info$variance.standardized <- sdnor(data,
                                            hvf.info$mean,
                                            sqrt(hvf.info$variance.expected),
                                            clip.max)
    #hvf.info$vst.variable <- order(hvf.info$variance.standardized) <2001
    topgenes <- hvf.info %>% top_n(nfeatures, variance.standardized)
    hvf.info$variable <- rownames(hvf.info) %in% rownames(topgenes)
    colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  }else{
    ###2. mean.var.plot method
    data <- object@assays[[assay]]@data
    mean.function = FastExpMean
    dispersion.function = FastLogVMR
    if (!inherits(x = mean.function, what = 'function')) {
      stop("'mean.function' must be a function")
    }
    if (!inherits(x = dispersion.function, what = 'function')) {
      stop("'dispersion.function' must be a function")
    }
    ##Calculate the variance to mean ratio (VMR) in non-logspace (return answer in log-space)
    feature.mean <- mean.function(data)
    feature.dispersion <- dispersion.function(data)
    names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = data)
    feature.dispersion[is.na(x = feature.dispersion)] <- 0
    feature.mean[is.na(x = feature.mean)] <- 0
    data.x.breaks <- switch(
      EXPR = binning.method,
      'equal_width' = num.bin,
      'equal_frequency' = c(
        -1,
        quantile(
          x = feature.mean[feature.mean > 0],
          probs = seq.int(from = 0, to = 1, length.out = num.bin)
        )
      ),
      stop("Unknown binning method: ", binning.method)
    )
    data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
    names(x = data.x.bin) <- names(x = feature.mean)
    mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
    sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
    feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
      sd.y[as.numeric(x = data.x.bin)]
    names(x = feature.dispersion.scaled) <- names(x = feature.mean)
    hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
    rownames(x = hvf.info) <- rownames(x = data)
    colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
    #cut off
    means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
    dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
    topfeatures <- rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    hvf.info$mvp.variable <-  means.use & dispersions.use
    if (selection.method == 'dispersion') {
      ##3. dispersion method
      hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), , drop = FALSE]
      hvf.info$disp.variable <- rownames(x = hvf.info) %in% head(x = rownames(x = hvf.info), n = nfeatures)
    }
  }
  return(hvf.info)
}

#' Find variable genes
#' @examples
#' data <- object@@assays[[assay]]@@counts
#' hvg.info <- FindVariableFeatures2(pbmc,selection.method = 'vst')
#' data <- object@@assays[[assay]]@@data
#' hvg.info <- FindVariableFeatures2(pbmc,selection.method = 'mean.var.plot')
#' @export
FindVariableFeatures2 <- function(data,
                                 selection.method = c('vst','mean.var.plot','dispersion'),
                                 loess.span = 0.3,
                                 clip.max = 'auto',
                                 num.bin = 20,
                                 binning.method = c('equal_width', 'equal_frequency'),
                                 nfeatures = 2000,
                                 mean.cutoff = c(0.1, 8),
                                 dispersion.cutoff = c(1, Inf)
) {
  selection.method <- match.arg(selection.method)
  binning.method <- match.arg(binning.method)
  if (selection.method == 'vst') {
    ###1. vst method.
    ##mean
    hvf.info <- data.frame(mean = Matrix::rowMeans(x = data))
    ##variance
    #hvf.info$variance <- rowVars(as.matrix(data))
    #as.matrix may cause memory problem when used for large data.
    hvf.info$variance <- apply(data, 1, var)
    #variance.expected
    hvf.info$variance.expected <- 0
    hvf.info$variance.standardized <- 0
    not.const <- hvf.info$variance > 0
    fit <- loess(
      formula = log10(x = variance) ~ log10(x = mean),
      data = hvf.info[not.const, ],
      span = loess.span
    )
    hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
    #variance.standardized
    sdnor <- function(x, mean.x, var.ex, clip.max){
      nor_x <- (x - mean.x)/var.ex
      #when value in var.ex equal to 0 ,return na. Then convert na to 0.
      nor_x[is.na(nor_x)] <- 0
      nor_x[nor_x > clip.max] <- clip.max
      sdx <- apply(nor_x, 1, var)
      return(sdx)
    }
    if (clip.max == 'auto') {
      clip.max <- sqrt(x = ncol(x = data))
    }
    hvf.info$variance.standardized <- sdnor(data,
                                            hvf.info$mean,
                                            sqrt(hvf.info$variance.expected),
                                            clip.max)
    #hvf.info$vst.variable <- order(hvf.info$variance.standardized) <2001
    topgenes <- hvf.info %>% top_n(nfeatures, variance.standardized)
    hvf.info$variable <- rownames(hvf.info) %in% rownames(topgenes)
    colnames(x = hvf.info) <- paste0('vst.', colnames(x = hvf.info))
  }else{
    ###2. mean.var.plot method
    mean.function = FastExpMean
    dispersion.function = FastLogVMR
    if (!inherits(x = mean.function, what = 'function')) {
      stop("'mean.function' must be a function")
    }
    if (!inherits(x = dispersion.function, what = 'function')) {
      stop("'dispersion.function' must be a function")
    }
    ##Calculate the variance to mean ratio (VMR) in non-logspace (return answer in log-space)
    feature.mean <- mean.function(data)
    feature.dispersion <- dispersion.function(data)
    names(x = feature.mean) <- names(x = feature.dispersion) <- rownames(x = data)
    feature.dispersion[is.na(x = feature.dispersion)] <- 0
    feature.mean[is.na(x = feature.mean)] <- 0
    data.x.breaks <- switch(
      EXPR = binning.method,
      'equal_width' = num.bin,
      'equal_frequency' = c(
        -1,
        quantile(
          x = feature.mean[feature.mean > 0],
          probs = seq.int(from = 0, to = 1, length.out = num.bin)
        )
      ),
      stop("Unknown binning method: ", binning.method)
    )
    data.x.bin <- cut(x = feature.mean, breaks = data.x.breaks)
    names(x = data.x.bin) <- names(x = feature.mean)
    mean.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = mean)
    sd.y <- tapply(X = feature.dispersion, INDEX = data.x.bin, FUN = sd)
    feature.dispersion.scaled <- (feature.dispersion - mean.y[as.numeric(x = data.x.bin)]) /
      sd.y[as.numeric(x = data.x.bin)]
    names(x = feature.dispersion.scaled) <- names(x = feature.mean)
    hvf.info <- data.frame(feature.mean, feature.dispersion, feature.dispersion.scaled)
    rownames(x = hvf.info) <- rownames(x = data)
    colnames(x = hvf.info) <- paste0('mvp.', c('mean', 'dispersion', 'dispersion.scaled'))
    #cut off
    means.use <- (hvf.info[, 1] > mean.cutoff[1]) & (hvf.info[, 1] < mean.cutoff[2])
    dispersions.use <- (hvf.info[, 3] > dispersion.cutoff[1]) & (hvf.info[, 3] < dispersion.cutoff[2])
    topfeatures <- rownames(x = hvf.info)[which(x = means.use & dispersions.use)]
    hvf.info$mvp.variable <-  means.use & dispersions.use
    if (selection.method == 'dispersion') {
      ##3. dispersion method
      hvf.info <- hvf.info[order(hvf.info$mvp.dispersion, decreasing = TRUE), , drop = FALSE]
      hvf.info$disp.variable <- rownames(x = hvf.info) %in% head(x = rownames(x = hvf.info), n = nfeatures)
    }
  }
  return(hvf.info)
}
