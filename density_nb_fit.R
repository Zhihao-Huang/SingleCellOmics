##var and mean
varmean_fun <- function (x, phi) {
  return(x + phi*x^2)
}
##no-zero fraction and mean
prob_zero_fun <- function (mu, phi) {
  if (phi == 0) {
    return(exp(-mu))
  }
  phi_1 = 1 / phi
  return((phi_1 / (mu + phi_1)) ** phi_1)
}

model_nb <- function(object, cluster, cut.off = 0.3, residual = .Machine$double.eps, 
                     assays = 'RNA', slot = 'counts') {
  message(cluster)
  cells <- names(Idents(object))[Idents(object) == cluster]
  subs <- subset(object, cells = cells)
  testdata <- subs@assays$RNA@counts
  Mean <- Matrix::rowMeans(testdata)
  Var <- apply(testdata, 1, var)
  zero_frac <- apply(testdata, 1, function(x) sum(x == 0)/length(x))
  x <- as.vector(Mean)
  y <- as.vector(Var)
  z <- as.vector(zero_frac)
  #plot(log10(x), log10(y))
  ##fit formula by python scipy.optimize.curve_fit
  phi <- fit_curve(x, y)
  ##var and mean
  pred <- x + phi*x^2
  plotdata <- data.frame(x = log10(x+residual), y = log10(y + residual), 
                         pred = log10(pred + residual),
                         z = z, hvg = rownames(subs) %in% subs@assays$RNA@var.features, 
                         mean = x, var = y,
                         Gene.name = rownames(testdata),phi = phi, stringsAsFactors = F)
  p1 <- ggplot(plotdata,aes(x,y)) + geom_point(size = 0.3) + geom_line(aes(x,pred), color = 'red', size = 0.5) + theme_bw() +
    xlab('log10(mean)') + ylab('log10(Var)')
  p2 <- ggplot(plotdata,aes(x,y,color = hvg)) + geom_point(size = 0.3) + geom_line(aes(x,pred), color = 'black', size = 0.5) + theme_bw() +
    xlab('log10(mean)') + ylab('log10(Var)')+ scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set1'))
  plotdata$prednb <- prob_zero_fun(x, phi)
  plotdata$preddeg <- plotdata$z > plotdata$prednb & plotdata$z < (1 - cut.off) 
  p <- ggplot(plotdata,aes(x,z,color = hvg)) + geom_point() + geom_line(aes(x,prednb), color = 'black', size = 1) + 
    geom_text(aes(label=ifelse(preddeg=='TRUE', Gene.name,"")), 
              size = 3, vjust=0.5, hjust= -0.2, position = position_dodge(0.5)) +
    theme_bw() + ggtitle(cluster) +
    xlab('log10(mean)') + ylab('Zero fraction') + scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set1'))
  message(paste(plotdata$Gene.name[plotdata$preddeg], collapse = ', '))
  return(plotdata)
}

library("reticulate")
library(ggplot2)
# Set the path to the Python executable file
#use_python("~/bin/python", required = T)
# Check the version of Python.
#py_config()
source_python('./optimize_curve_fit.py')
#Load data. https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
pbmcmat <- Read10X(data.dir = './filtered_gene_bc_matrices/hg19/')
dim(pbmcmat)
#[1] 32738  2700
pbmc <- CreateSeuratObject(counts = pbmcmat, min.cells = 3)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
dim(pbmc)
#[1] 13714  2638
pbmc <- NormalizeData(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'dispersion', nfeatures = 2000)
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
pbmc <- RunPCA(pbmc)
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, dims = 1:10, resolution = 0.5)
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, label = T)
FeaturePlot(pbmc,pt.size = 0.1, features = c("MS4A1", "GNLY", "CCR7", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", 
                                             "CD8A", 'IL7R','CD4','S100A4'))
new.cluster.ids <- c("Naive CD4 T","Memory CD4 T",  "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc)
pbmc@meta.data$Annotation <- Idents(pbmc)

#fit model
datalist <- list()
for (i in unique(Idents(pbmc))) {
  data <- model_nb(pbmc, i, cut.off = 0, residual = 0.00001)
  datalist[[i]] <- data
}

##in groups
plist <- list()
for (i in unique(Idents(pbmc))) {
  p <- ggplot(datalist[[i]],aes(x,z,color = hvg)) + geom_point(size = 0.1) + 
    geom_line(aes(x,prednb), color = 'black', size = 1) + 
    theme_bw() + theme(legend.position = "none") + ggtitle(i) +
    xlab('log10(mean)') + ylab('Zero fraction') + scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set1'))
  plist[[i]] <- p
}
do.call(gridExtra::grid.arrange, plist)
#density
plist <- list()
for (i in unique(Idents(pbmc))) {
  data <- datalist[[i]]
  p <- ggplot(data[data$hvg == 'TRUE',],aes(x,z)) + 
    geom_density_2d(size = 0.3)+
   # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    geom_line(aes(x,prednb), color = 'black', size = 1) + 
    theme_bw() + theme(legend.position = "none") + ggtitle(i) +
    xlab('log10(mean)') + ylab('Zero fraction')
  plist[[i]] <- p
}
do.call(gridExtra::grid.arrange, plist)

plist <- list()
for (method in c('vst','mean.var.plot','dispersion')) {
  pbmc <- FindVariableFeatures(pbmc, selection.method = method, nfeatures = 2000)
  Idents(pbmc) <- 1
  plotdata <- model_nb(pbmc, 1, residual = 0.00001)
  Idents(pbmc) <- pbmc$Annotation
  p1 <- ggplot(plotdata,aes(x,z,color = hvg)) + geom_point(size = 0.1) + 
    geom_line(aes(x,prednb), color = 'black', size = 1) + 
    theme_bw() + theme(legend.position = "none") + ggtitle(method) +
    xlab('log10(mean)') + ylab('Zero fraction') + scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set1'))
  p2 <- ggplot(plotdata[plotdata$hvg == 'TRUE',],aes(x,z)) + 
    geom_density_2d(size = 0.3)+
    # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    geom_line(aes(x,prednb), color = 'black', size = 1) + 
    theme_bw() + theme(legend.position = "none") +
    xlab('log10(mean)') + ylab('Zero fraction')
  plist[[method]] <- p1 + p2
}

plist[[2]]
#all genes density
ggplot(plotdata,aes(x,z)) + 
  geom_density_2d(size = 0.3)+
  # stat_density_2d(aes(fill = ..level..), geom = "polygon") +
  geom_line(aes(x,prednb), color = 'black', size = 1) + 
  theme_bw() +  ggtitle('All genes') +
  xlab('log10(mean)') + ylab('Zero fraction')
ggplot(plotdata,aes(x,z,color = hvg)) + geom_point(size = 0.1) + 
  geom_line(aes(x,prednb), color = 'black', size = 1) + 
  theme_bw() + ggtitle('All genes') +
  xlab('log10(mean)') + ylab('Zero fraction') + scale_color_manual(values = RColorBrewer::brewer.pal(3,'Set1'))

#HVG enriched in over-dispersion region.
datamerge_list <- list()
for (method in c('vst','mean.var.plot','dispersion')) {
  pbmc <- FindVariableFeatures(pbmc, selection.method = method, nfeatures = 2000)
  datalist <- list()
  for (i in unique(Idents(pbmc))) {
    data <- model_nb(pbmc, i, cut.off = 0, residual = 0.00001)
    datalist[[i]] <- data
  }
  plotdatalist <- list()
  for (i in unique(Idents(pbmc))) {
    data <- datalist[[i]]
    over_index <- data$y > data$pred 
    datao <- data[over_index,]
    datao <- data.frame(perc = table(datao$hvg)[2]/table(datao$hvg)[1], cell = i)
    plotdatalist[[i]] <- datao
  }
  plotdata_merge <- do.call(rbind, plotdatalist)
  plotdata_merge$method <- method
  datamerge_list[[method]] <- plotdata_merge
}
plotdata_merge <- do.call(rbind,datamerge_list)
plotdata_merge$method <- factor(plotdata_merge$method, levels =  c('dispersion','vst','mean.var.plot'))
ggplot(plotdata_merge, aes(cell, perc, fill = method)) + geom_bar(width = 0.5, stat = 'identity',position = 'dodge') + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  geom_hline(yintercept=2000/nrow(data), linetype="dashed") +
  xlab('') + ylab('Percentage of HVG in low-detection_rate zone')


##gene bin
#group by mean expression
feature.mean <- apply(pbmc@assays$RNA@data,1,mean)
data.x.bin <- cut(x = feature.mean, breaks = 9)
names(x = data.x.bin) <- rownames(pbmc)
pbmc <- FindVariableFeatures(pbmc, selection.method = 'dispersion', nfeatures = 2000)
meta <- pbmc@assays$RNA@meta.features
meta$bin <- as.numeric(data.x.bin)
meta$dsp.variable <- rownames(meta) %in% pbmc@assays$RNA@var.features
# measure value of dispersion
id.var <-c('mvp.dispersion','vst.variance.standardized', 'mvp.dispersion.scaled')
meta2 <- meta[,c(id.var,'bin')]
colnames(meta2) <- c('dispersion','vst','mean.var.plot','Gene groups')
data_bin <- melt(meta2, id.vars ='Gene groups')
data_bin$bin <- as.character(data_bin$`Gene groups`)
ggplot(data_bin, aes(bin, value, fill = variable)) + geom_boxplot(outlier.size = 0.1) + 
  theme_bw()+ xlab('Gene groups') + ylab('Value of dispersion')
# measure number of HVG per bin
sum.num <- data.frame(vst = tapply(X = meta$vst.variable, INDEX = meta$bin, FUN = sum),
                      dispersion = tapply(X = meta$dsp.variable, INDEX = meta$bin, FUN = sum),
                      mean.var.plot= tapply(X = meta$mvp.variable, INDEX = meta$bin, FUN = sum))
sum.num$bin <- rownames(sum.num)
data_bin <- melt(sum.num, id.vars = 'bin')
data_bin$variable <- factor(data_bin$variable, levels =  c('dispersion','vst','mean.var.plot'))
ggplot(data_bin, aes(bin, value, fill = variable)) + geom_bar(stat = 'identity', position = 'dodge') + 
  theme_bw()+ xlab('Gene groups') + ylab('Number of HVG')
