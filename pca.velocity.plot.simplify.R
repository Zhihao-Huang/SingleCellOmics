pca.velocity.plot.simplify <- function(resultlist, pc.multipliers = c(1,-1), 
                                       arrow.scale = 1, nPcs = 2, show.grid.flow = F, 
                                       grid.n = 20, grid.sd = NULL, pcount = 1, arrow.lwd = 1,
                                       size.norm = FALSE, cell.colors = NULL, plot.cols=1,
                                       plot.grid.points = FALSE, fixed.arrow.length = FALSE, 
                                       max.grid.arrow.length = NULL,min.arrow.size = NULL,
                                       min.grid.cell.mass = 1) {
  x0 <- resultlist$current
  x1 <- resultlist$projected
  message("Step1. log ... ")
  x0.log <- log2(x0 + pcount)
  x1.log <- log2(x1 + pcount)
  message("Step2. pca ... ")
  cent <- rowMeans(x0.log)
  epc <- pcaMethods::pca(t(x0.log - cent), center = F, nPcs = nPcs)
  message("Step3. pc multipliers ... ")
  epc@loadings <- t(t(epc@loadings) * pc.multipliers)
  epc@scores <- scale(epc@completeObs, scale = F, center = T) %*% epc@loadings
  message('Step4. projection of the extrapolated state on the same eigenvectors...  ')
  x1.scores <- t(x1.log - cent) %*% epc@loadings
  message("Step5. delta norm , scale arrow size... ")
  delta.pcs <- as.matrix(x1.scores - epc@scores)
  delta.pcs <- delta.pcs * arrow.scale
  x1.scores <- epc@scores + delta.pcs
  message("Step6. plot arrow...")
  par(mfrow = c(ceiling((nPcs - 1)/plot.cols), plot.cols), 
      mar = c(3.5, 3.5, 2.5, 1.5), mgp = c(2, 0.65, 0), 
      cex = 0.85)
  vinfo <- lapply(1:(nPcs - 1), function(i) {
    pos <- epc@scores[, c((i - 1) + 1, (i - 1) + 2)]
    ppos <- pos + delta.pcs[, c((i - 1) + 1, (i - 1) + 2)]
    plot(pos, bg = cell.colors[rownames(pos)], pch = 21, 
         col = ac(1, alpha = 0.3), lwd = 0.5, 
         xlab = paste("PC", (i - 1) + 1), ylab = paste("PC", (i - 1) + 2), 
         axes = T, 
         main = paste("PC", (i - 1) + 1, " vs. PC", (i - 1) + 2, sep = ""))
    box()
    if (show.grid.flow) {
      ars <- data.frame(pos[, 1], pos[, 2], ppos[, 1], ppos[, 2])
      colnames(ars) <- c("x0", "y0", "x1", "y1")
      arsd <- data.frame(xd = ars$x1 - ars$x0, yd = ars$y1 - ars$y0)
      rownames(ars) <- rownames(arsd) <- rownames(pos)
      rx <- range(c(range(ars$x0), range(ars$x1)))
      ry <- range(c(range(ars$y0), range(ars$y1)))
      gx <- seq(rx[1], rx[2], length.out = grid.n)
      gy <- seq(ry[1], ry[2], length.out = grid.n)
      if (is.null(grid.sd)) {
        grid.sd <- sqrt((gx[2] - gx[1])^2 + (gy[2] - gy[1])^2)/2
        cat("grid.sd=", grid.sd, " ")
      }
      if (is.null(min.arrow.size)) {
        min.arrow.size <- sqrt((gx[2] - gx[1])^2 + (gy[2] - gy[1])^2) * 0.01
        cat("min.arrow.size=", min.arrow.size, " ")
      }
      if (is.null(max.grid.arrow.length)) {
        max.grid.arrow.length <- sqrt(sum((par("pin")/c(length(gx), length(gy)))^2)) * 0.25
        cat("max.grid.arrow.length=", max.grid.arrow.length, " ")
      }
      garrows <- do.call(rbind, lapply(gx, function(x) {
        cd <- sqrt(outer(pos[, 2], -gy, "+")^2 + (x - pos[, 1])^2)
        cw <- dnorm(cd, sd = grid.sd)
        gw <- Matrix::colSums(cw)
        cws <- pmax(1, Matrix::colSums(cw))
        gxd <- Matrix::colSums(cw * arsd$xd)/cws
        gyd <- Matrix::colSums(cw * arsd$yd)/cws
        al <- sqrt(gxd^2 + gyd^2)
        vg <- gw >= min.grid.cell.mass & al >= min.arrow.size
        cbind(rep(x, sum(vg)), gy[vg], x + gxd[vg], gy[vg] + gyd[vg])
      }))
      colnames(garrows) <- c("x0", "y0", "x1", "y1")
      if (fixed.arrow.length) {
        suppressWarnings(arrows(garrows[, 1], garrows[,2], garrows[, 3], garrows[, 4], length = 0.05, 
                                lwd = arrow.lwd))
      }
      else {
        alen <- pmin(max.grid.arrow.length, sqrt(((garrows[,3] - garrows[, 1]) * par("pin")[1]/diff(par("usr")[c(1,2)]))^2 + ((garrows[, 4] - garrows[, 2]) * par("pin")[2]/diff(par("usr")[c(3, 4)]))^2))
        suppressWarnings(lapply(1:nrow(garrows), function(i) arrows(garrows[i,1], garrows[i, 2], garrows[i, 3], garrows[i, 4], length = alen[i], lwd = arrow.lwd)))
      }
      if (plot.grid.points) 
        points(rep(gx, each = length(gy)), rep(gy, length(gx)), 
               pch = ".", cex = 0.1, col = ac(1, alpha = 0.4))
    }
    else {
      grid()
      suppressWarnings(arrows(pos[, 1], pos[, 2], ppos[, 1], ppos[, 2], length = 0.05, lwd = arrow.lwd))
    }
  })
}
ac <- function (x, alpha = 1, ...) 
{
  y <- adjustcolor(x, alpha.f = alpha, ...)
  names(y) <- names(x)
  return(y)
}
