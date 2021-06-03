source('velocity.estimates.simplify.R')
source('pca.velocity.plot.simplify.R')
#load data
emat <- read.table('emat.txt')
nmat <- read.table('nmat.txt')
comdf <- read.table('cell_cluster.txt', header = T)
cluster.label <- factor(comdf$cluster)
names(cluster.label) <- comdf$cellID
cell.colors = MUDAN:::fac2col(cluster.label)
## RNA velocity model without pooling
rvel.cd.unpooled <- velocity.estimates.simplify(emat, nmat,
                                                     fit.quantile = 0.05,
                                                     min.nmat.emat.correlation = 0.2, 
                                                     min.nmat.emat.slope = 0.2)
## plot onto PCA
#unpooled
pca.velocity.plot.simplify(rvel.cd.unpooled,
                  nPcs=2,
                  plot.cols=1,
                  cell.colors=cell.colors,
                  pc.multipliers=c(1,-1) ## adjust as needed to orient pcs
)
pca.velocity.plot.simplify(rvel.cd.unpooled,
                  nPcs=2,
                  plot.cols=1,
                  cell.colors=cell.colors,
                  pc.multipliers=c(1,-1), ## adjust as needed to orient pcs
                  show.grid.flow = TRUE, 
                  grid.n=20 ## adjust as needed
)

###plot cell-cycle phase
library(Seurat)
library(ggplot2)
rownames(comdf) <- comdf$cellID
seus4 <- CreateSeuratObject(counts = emat, meta.data = comdf)
seus4 <- CellCycleScoring(object = seus4, g2m.features = cc.genes$g2m.genes, 
                                  s.features = cc.genes$s.genes)
#pca
x0 <- rvel.cd.unpooled$current
x1 <- rvel.cd.unpooled$projected
message("Step1. log ... ")
x0.log <- log2(x0 + 1)
x1.log <- log2(x1 + 1)
message("Step2. pca ... ")
cent <- rowMeans(x0.log)
epc <- pcaMethods::pca(t(x0.log - cent), center = F, nPcs = 2)
message("Step3. pc multipliers ... ")
epc@loadings <- t(t(epc@loadings) * c(1,-1))
epc@scores <- scale(epc@completeObs, scale = F, center = T) %*% epc@loadings
head(epc@scores)
ccdf <- as.data.frame(epc@scores)
ccdf <- cbind(ccdf, seus4@meta.data[,c('S.Score','G2M.Score','Phase')])
ggplot(ccdf, aes(PC1, PC2, color = Phase)) + geom_point() + theme_bw() 
ggplot(ccdf, aes(PC1, PC2, color = S.Score)) + geom_point() + theme_bw() +
  scale_color_gradientn(colours = c('grey','yellow','red'))
ggplot(ccdf, aes(PC1, PC2, color = G2M.Score)) + geom_point() + theme_bw() +
  scale_color_gradientn(colours = c('grey','yellow','red'))
