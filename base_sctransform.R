#' basic method of sctransform
#' 
#' @examples 
#' datalist <- sctransform_base(pbmc@@assays$RNA@@counts, smooth_th = 100, min.cell = 5)
sctransform_base <- function(raw.count, return.data = c('correct_mat', 'list')[2], min_cells = 5,
                             smooth_th = 10) {
  #####sctransform by base function
  ####Step1. Fit independent regression models per gene.
  #a.Get log10 total cell umi-count.
  umi <- raw.count
  log10_total_umi <- log10(apply(umi, 2, sum)) 
  genes_cell_count <- apply(umi >= 0.01, 1, sum)
  min_cells <- min_cells
  genes <- rownames(umi)[genes_cell_count >= min_cells]
  umi <- umi[genes, ]
  #b. Fit regression per gene. This step is time consuming.
  #The default method in Seurat::SCTransform is  quasipoisson. 
  message('Step1. Fit regression per gene.')
  fitlist <- pbapply::pblapply(1:nrow(umi), function(x) {
    one_gene_umi <- umi[x,]
    df <- data.frame(x = one_gene_umi, y = log10_total_umi, stringsAsFactors = F)
    #fit quasi-poisson
    fit <- glm(x~y, family = quasipoisson, data = df)
    #estimate theta of the Negative Binomial Distribution
    theta <- as.numeric(x = MASS::theta.ml(y = one_gene_umi, mu = fit$fitted))
    c(theta, fit$coefficients)
  }) 
  fitmat <- do.call(rbind, fitlist)
  colnames(fitmat) <- c('theta','coef1','coef2')
  rownames(fitmat) <- rownames(umi)
  min_theta <- 1e-7
  fitmat[, 'theta'] <- pmax(fitmat[, 'theta'], min_theta)
  ####Step2. Exploit the relationship of model parameter values and gene mean to learn global trends.
  # variance of NB is mu * (1 + mu / theta)
  # (1 + mu / theta) is what we call overdispersion factor here
  # a. Use overdispersion factor instead of theta.
  message('Step2. Exploit the relationship of model parameter values and gene mean to learn global trends.')
  genes_log_gmean <- log10(exp(apply(log(umi+1),1,mean)) - 1)
  dispersion_par <-  log10(1 + 10^genes_log_gmean / fitmat[, 'theta'])
  model_pars <- fitmat[, colnames(fitmat) != 'theta']
  model_pars <- cbind(dispersion_par, model_pars)
  # b. look for outliers in the parameters
  # outliers are those that do not fit the overall relationship with the mean at all
  is_outlier <- function(y, x, th = 10) {
    bin.width <- (max(x) - min(x)) * bw.SJ(x) / 2
    eps <- .Machine$double.eps * 10
    breaks1 <- seq(from = min(x) - eps, to = max(x) + bin.width, by = bin.width)
    breaks2 <- seq(from = min(x) - eps - bin.width/2, to = max(x) + bin.width, by = bin.width)
    score1 <- robust_scale_binned(y, x, breaks1)
    score2 <- robust_scale_binned(y, x, breaks2)
    return(pmin(abs(score1), abs(score2)) > th)
  }
  robust_scale_binned <- function(y, x, breaks) {
    bins <- cut(x = x, breaks = breaks, ordered_result = TRUE)
    tmp <- aggregate(x = y, by = list(bin=bins), FUN = robust_scale)
    score <- rep(0, length(x))
    o <- order(bins)
    if (inherits(x = tmp$x, what = 'list')) {
      score[o] <- unlist(tmp$x)
    } else {
      score[o] <- as.numeric(t(tmp$x))
    }
    return(score)
  }
  robust_scale <- function(x) {
    return((x - median(x)) / (mad(x) + .Machine$double.eps))
  }
  outliers <- apply(model_pars, 2, function(y) is_outlier(y, genes_log_gmean, th = smooth_th))
  outliers <- apply(outliers, 1, any)
  genes <- rownames(model_pars)
  model_pars <- model_pars[!outliers, ]
  genes_log_gmean_filtered <- genes_log_gmean[!outliers]
  umi_filtered <- umi[!outliers,]
  
  # c. Select bandwidth to be used for smoothing
  bw <- bw.SJ(genes_log_gmean_filtered) * bw_adjust
  # d. For parameter predictions
  x_points <- pmax(genes_log_gmean, min(genes_log_gmean_filtered))
  x_points <- pmin(x_points, max(genes_log_gmean_filtered))
  # e. Take results from step 1 and fit/predict parameters to all genes
  o <- order(x_points)
  model_pars_fit <- matrix(NA_real_, length(genes), ncol(model_pars),
                           dimnames = list(genes, colnames(model_pars)))
  # f. fit / regularize dispersion parameter
  model_pars_fit[o, 'dispersion_par'] <- ksmooth(x = genes_log_gmean_filtered, y = model_pars[, 'dispersion_par'],
                                                 x.points = x_points, bandwidth = bw, kernel='normal')$y
  # g. global fit / regularization for all coefficients
  for (i in 2:ncol(model_pars)) {
    model_pars_fit[o, i] <- ksmooth(x = genes_log_gmean_filtered, y = model_pars[, i],
                                    x.points = x_points, bandwidth = bw, kernel='normal')$y
  }
  # h. back-transform dispersion parameter to theta
  theta <- 10^genes_log_gmean / (10^model_pars_fit[, 'dispersion_par'] - 1)
  model_pars_fit <- model_pars_fit[, colnames(model_pars_fit) != 'dispersion_par']
  model_pars_fit <- cbind(theta, model_pars_fit)
  ####Step3. Use the regularized regression parameters to define an affine function that
  #transforms UMI counts into Pearson residuals.
  message('Step3. Use the regularized regression parameters to define an affine function that transforms UMI counts into Pearson residuals.')
  coefs <- model_pars_fit[, -1, drop=FALSE]
  theta <- model_pars_fit[, 1]
  # get pearson residuals
  regressor_data_orig <- cbind(rep(1,length(log10_total_umi)), log10_total_umi)
  mu <- exp(tcrossprod(coefs, regressor_data_orig))
  variance <- mu + mu^2 / theta
  pearson_residual <- (umi - mu) / sqrt(variance)
  # generate output
  regressor_data_median <- regressor_data_orig
  regressor_data_median[,'log10_total_umi'] <- apply(regressor_data_median[, 'log10_total_umi', drop=FALSE], 2, function(x) rep(median(x), length(x)))
  mu <- exp(tcrossprod(coefs, regressor_data_median))
  variance <- mu + mu^2 / theta
  y.res <- mu + pearson_residual * sqrt(variance)
  y.res <- round(y.res, 0)
  y.res[y.res < 0] <- 0
  correct_mat <- as(y.res, Class = 'dgCMatrix')
  umi_corrected_log1p <- as(log(correct_mat + 1), Class = 'dgCMatrix')
  #top var gene
  pr <- pearson_residual
  clip.range <- c(-sqrt(ncol(raw.count)), sqrt(ncol(raw.count)))
  pr[pr < clip.range[1]] <- clip.range[1]
  pr[pr > clip.range[2]] <- clip.range[2]
  gene_var <- apply(pr, 1, var)
  gene_var_sorted <- gene_var[order(gene_var, decreasing = T)]
  # clip the residuals for scale.data
  scale.data <- pearson_residual
  # default: sctransform: sqrt(ncol(umi); Seurat: sqrt(ncol(umi)/30
  clip.range <- c(-sqrt(ncol(raw.count) / 30), sqrt(ncol(raw.count) / 30))
  scale.data[scale.data < clip.range[1]] <- clip.range[1]
  scale.data[scale.data > clip.range[2]] <- clip.range[2]
  
  # 2nd regression
  scale.data <- Seurat::ScaleData(
    scale.data,
    features = NULL,
    vars.to.regress = NULL,
    latent.data = NULL,
    model.use = 'linear',
    use.umi = FALSE,
    do.scale = F,
    do.center = T,
    scale.max = Inf,
    block.size = 750,
    min.cells.to.block = 3000,
    verbose = T
  )
  datalist <- list('count' = correct_mat,
                   'data' = umi_corrected_log1p,
                   'model_pars_fit' = model_pars_fit,
                   'model_pars_outliers' = outliers,
                   'gene_var' = gene_var_sorted,
                   'top.features' = names(gene_var_sorted),
                   'scale.data' = scale.data)
  data <- switch(return.data,
                 'correct_mat' = correct_mat,
                 'list' = datalist,
                 stop('data ', return.data, ' unknown - only correct_mat and list supported at the moment.')
  )
  return(data)
}

pbmc <- readRDS('./pbmc.rds')
datalist <- sctransform_base(pbmc@assays$RNA@counts, smooth_th = 10, min_cells = 5)
library(Seurat)
assay.out <- CreateAssayObject(counts = datalist$count)
assay.out <- SetAssayData(
  object = assay.out,
  slot = 'data',
  new.data = datalist$data
)
assay.out <- SetAssayData(
  object = assay.out,
  slot = 'scale.data',
  new.data = datalist$scale.data
)
#assay.out[['sct.variable']] <- rownames(x = assay.out[[]]) %in% datalist$top.features[1:3000]
pbmc[['SCTb']] <- assay.out
pbmc[['SCTb']]@var.features <- datalist$top.features[1:3000]
DefaultAssay(pbmc) <- 'SCTb'
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:30)
pbmc <- FindClusters(pbmc, dims = 1:30, resolution = 1.8)
pbmc <- RunUMAP(pbmc, dims = 1:30)
DimPlot(pbmc, label = T)
DimPlot(pbmc, group.by = 'Annotation')