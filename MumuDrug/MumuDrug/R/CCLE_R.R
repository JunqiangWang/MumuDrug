


#' @export
StochDownSamp <- function(count.mat, depth.vect) {

  require(DirichletReg)
  require(Matrix)
  require(optparse)


  # create empty down sampled matrix
  ds.mat <- matrix(0L, nrow = nrow(count.mat), ncol = ncol(count.mat))
  rownames(ds.mat) <- rownames(count.mat); colnames(ds.mat) <- colnames(count.mat)
  n.cells <- ncol(count.mat)
  # make downsampled cells
  for (i in 1:n.cells) {
    p.vect <- count.mat[,i]; p.vect <- p.vect / sum(p.vect)
    ds.vect <- rmultinom(1, depth.vect[i], prob = p.vect)
    ds.mat[,i] <- ds.vect
  }
  # return matrix
  return(ds.mat)
}

#




#' @export
CPMTransform <- function(dat.mat, l2 = FALSE, pseudo = FALSE) {
  if (pseudo) {
    dat.mat <- dat.mat + 1
  }
  cpm.mat <- t(t(dat.mat) / (colSums(dat.mat) / 1e6))
  if (l2) {
    cpm.mat <- log2(cpm.mat + 1)
  }
  return(cpm.mat)
}


#

#' @export
SignatureByDownSample<-function(sc.mat, ccle.mat){


  shared.genes <- intersect(rownames(sc.mat), rownames(ccle.mat))
  ccle.mat <- ccle.mat[shared.genes,]
  sc.mat <- sc.mat[shared.genes,]
  ds.depth.vec <- rep(mean(colSums(sc.mat)), ncol(ccle.mat))


  ccle.ds <- StochDownSamp(ccle.mat, ds.depth.vec)

  ccle.ss <- ccle.ds[which(rowSums(ccle.ds) > 0),]

  ccle.ss <- CPMTransform(ccle.ss)

  # intersect genes
  shared.genes <- intersect(rownames(sc.mat), rownames(ccle.ss))

  sc.mat <- sc.mat[shared.genes,]
  ccle.ss <- ccle.ss[shared.genes,]

  # normalize
  ccle.mean <- rowMeans(ccle.ss)
  ccle.sd <- apply(ccle.ss, 1, sd)
  sc.mat <- apply(sc.mat, 2, function(x) { (x - ccle.mean) / ccle.sd } )

  return(sc.mat)

}



