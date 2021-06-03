velocity.estimates.simplify <- function(emat,nmat,deltaT=1,steady.state.cells=colnames(emat),
                                        mult=1e3,min.nmat.smat.correlation=0.05,
                                        min.nmat.emat.correlation=0.05, min.nmat.emat.slope=0.05, 
                                        zero.offset=FALSE,deltaT2=1, fit.quantile=NULL, 
                                        diagonal.quantiles=FALSE, kCells = 1,
                                        do.par=TRUE,  n.cores=4) {
  resl <- list();
  # bring matrices to the same gene set (just in case)
  vg <- intersect(rownames(emat),rownames(nmat));
  emat <- emat[vg,]; nmat <- nmat[vg,]
  pcount <- 1;
  # size estimates
  emat.size <- Matrix::colSums(emat)
  nmat.size <- Matrix::colSums(nmat)
  emat.cs <- emat.size[colnames(emat)]/mult;
  nmat.cs <- nmat.size[colnames(nmat)]/mult;
  emat.log.norm <- log(as.matrix(t(t(emat)/emat.cs))+pcount);
  conv.emat <- emat
  conv.nmat <- nmat
  conv.emat.cs <- emat.cs
  conv.nmat.cs <- nmat.cs
  # size-normalized counts
  conv.emat.norm <- t(t(conv.emat)/conv.emat.cs)
  conv.nmat.norm <- t(t(conv.nmat)/conv.nmat.cs)
  # size-normalized counts
  emat.norm <- t(t(emat)/emat.cs)
  nmat.norm <- t(t(nmat)/nmat.cs)
  resl$conv.nmat.norm <- conv.nmat.norm;
  resl$conv.emat.norm <- conv.emat.norm;
  cat("fitting gamma coefficients ... ")
  ko <- data.frame(do.call(rbind,parallel::mclapply(sn(rownames(conv.emat.norm)),function(gn) {
    df <- data.frame(n=(conv.nmat.norm[gn,steady.state.cells]),
                     e=(conv.emat.norm[gn,steady.state.cells]))
    o <- 0;
    #########################################
    ###equivalent to: 
    ###zi <- emat < 1; 
    ###proportion of cells that spliced (exonic) count < 1.
    ###o <-  sum(df$n[zi])/(sum(zi)+1)
    #########################################
    zi <- df$e<1/conv.emat.cs[steady.state.cells]; 
    if(any(zi)) { o <- sum(df$n[zi])/(sum(zi)+1)} 
    df$o <- o;
    eq <- quantile(df$e,p=c(fit.quantile,1-fit.quantile))
    pw <- as.numeric(df$e>=eq[2] | df$e<=eq[1])
    d <- lm(n~e,data=df,weights=pw)
    # note: re-estimating offset here
    return(c(o=as.numeric(coef(d)[1]),g=as.numeric(coef(d)[2]),r=cor(df$e,df$n,method='spearman')))
  },mc.cores=n.cores,mc.preschedule=T)))
  ko <- na.omit(ko)
  cat("done. successful fit for",nrow(ko),"genes\n")
  full.ko <- ko;
  vi <- ko$r>min.nmat.emat.correlation
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low nmat-emat correlation\n")
  ko <- ko[vi,]
  vi <- ko$g>min.nmat.emat.slope
  if(!all(vi)) cat("filtered out",sum(!vi),"out of",length(vi),"genes due to low nmat-emat slope\n")
  ko <- ko[vi,]
  gamma <- ko$g; 
  offset <- ko$o; 
  names(gamma) <- names(offset) <- rownames(ko);
  cat("calculating RNA velocity shift ... ")
  deltaE <- t.get.projected.delta(conv.emat.norm,conv.nmat.norm,gamma,offset=offset,delta=deltaT)
  resl$gamma <- gamma;
  cat("done\n")
  cat("calculating extrapolated cell state ... ")
  # reduced cell normalization (only genes for which momentum was estimated)
  emat.norm <- emat[rownames(emat) %in% rownames(deltaE),]
  emat.sz <- emat.cs;
  emat.norm <- t(t(emat.norm)/(emat.sz));
  ######
  ### calculates the difference in the number of counts based on the library size, renormalizes
  ### note: also introduces centripetal velocity
  ######
  emn <- t.get.projected.cell2(emat.norm,emat.sz,as.matrix(deltaE),mult = mult,delta=deltaT2);
  cat("done\n")
  full.ko$valid <- rownames(full.ko) %in% rownames(ko)
  resl <- c(resl,list(projected=emn,current=emat.norm,deltaE=deltaE,deltaT=deltaT,ko=full.ko,mult=mult,kCells=kCells));
  return(resl)
} 

# get gene strand information and classify concordance
sn <- function(x) { names(x) <- x; x}
# estimate projected delta given x'=(y-o) - gamma*x solution
# em - normalized expression matrix
# nm - normalized nascent matrix
# gamma - inferred degradation coefficients
# o - inferred offset (assumed to be zero by default)
# delta - time to project forward
t.get.projected.delta <- function(em,nm,gamma,offset=rep(0,length(gamma)),delta=0.5) {
  # adjust rownames
  gn <- intersect(names(gamma),rownames(em));
  if(is.null(names(offset))) { names(offset) <- names(gamma); }
  em <- em[gn,]; nm <- nm[gn,]; gamma <- gamma[gn]; offset <- offset[gn];
  # time effect constant
  egt <- exp(-gamma*delta);
  y <- nm-offset; y[y<0] <- 0; # zero out entries with a negative n levels after offset adjustment
  em*egt + (1-egt)*y/gamma  - em
}
# calculates the difference in the number of counts based on the library size, renormalizes
# note: also introduces centripetal velocity
t.get.projected.cell2 <- function(em,cellSize,deltae,mult=1e3,delta=1) {
  rz <- matrix(0,nrow=nrow(em),ncol=ncol(em)); colnames(rz) <- colnames(em); rownames(rz) <- rownames(em)
  gn <- intersect(rownames(deltae),rownames(rz))
  rz[match(gn,rownames(rz)),colnames(deltae)] <- deltae[gn,]; 
  # translate fpm delta into the number of molecules based on the current cell size
  rz <- t(t(rz)*cellSize)
  emm <- t(t(em)*cellSize)
  emn <- emm + rz*delta; 
  emn[emn<0] <- 0;
  newCellSize <- (cellSize+Matrix::colSums(emn-emm)/mult)
  emn <- t(t(emn)/newCellSize)
  
  #emn <- t(t(emn)/Matrix::colSums(emn)*Matrix::colSums(em))
  emn
}