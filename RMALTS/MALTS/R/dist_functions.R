# 
mahal_dist_wv <- function(X, z, w) {
  if (length(z) != nrow(X)) {
    stop("Length of z does not match the number of rows in X")
  }

  if (ncol(X) != length(w)) {
    stop("Number of features in X does not match then number of weights w")
  }
  n_1 <- sum(z)
  n_0 <- length(z) - n_1
  dist_mat <- matrix(NA, n_1, n_0)

  X_0 <- X[z == 0, , drop = FALSE]
  X_1 <- X[z == 1, , drop = FALSE]
  for (i in 1:n_1) {
    dist_mat[i, ] <- mahalanobis(X_0, X_1[i, ], cov = diag(w^2), inverted = TRUE)
  }

  return(dist_mat)
}


cross_dist <- function(X,Y,w) {
  Xs <- reweight_mat(X, w)
  Ys <- reweight_mat(Y,w)

  if(missing(Y)) Ys <- NULL
  d_mat <- proxy::dist(Xs,Ys)^2
  return(matrix(d_mat, nrow = nrow(d_mat), ncol = ncol(d_mat)))
}


cross_dist_rfast <- function(X,Y,w) {
  Xs <- reweight_mat(X, w)
  Ys <- reweight_mat(Y,w)

  if(missing(Y)) Ys <- NULL
  d_mat <- Rfast::dista(Xs,Ys, square = TRUE)
  return(d_mat)
}

cross_dist_rfast2 <- function(X,Y,w) {
  Xs <- reweight_mat(X, w)
  Ys <- reweight_mat(Y,w)

  if(missing(Y)) Ys <- NULL
  d_mat <- collapse::setop(
    Rfast::rowsums(Xs^2) - 2*tcrossprod(Xs, Ys), "+", Rfast::rowsums(Ys^2), 
		 rowwise = TRUE)
  return(d_mat)
}


reweight_mat <- function(X, w) {
  return(sweep(X, 2L, w, `*`))
}


