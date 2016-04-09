#' Shrunken covariance matrix
#'
#' Calculates a shrunken covariance matrix to be passed as the demoninator
#' in the LDA formula. The shrunken matrix is in the subspace of the components
#' that capture 95\% of the variation of the dataset.
#'
#' @param R PCA scores
#' @param samples_count number of samples to explain 95\% of the variance
#' @param r regularisation factor to smooth covariance matrix in order to
#'	reduce noise. 1 is no regularisatoin.
#' 
#' @return covariance matrix
#' @export

shrink_mat <- function(R, samples_count, r){

    Dd_pre <- t(R) %*% R / samples_count
    Dd <- diag(diag(Dd_pre))
    sigma <- mean(diag(Dd))
    shrunkMat <- r * Dd + sigma * (1 - r) * diag(dim(R)[2])
    return(shrunkMat)

}
