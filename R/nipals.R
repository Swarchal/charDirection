#' Principal component analysis
#' 
#' Calculates the principal components accorindg to the algorithm by Wold
#'
#' @param X matrix or numeric dataframe.
#' @param a integer, number of principal components to calculate.
#' @param it integer, maximum number of iterations if error tolerance is not
#'	reached.
#' @param tolerance float, precision tolerance.
#' 
#' @return list, T: PCA scores. P: PCA loadings. pcvar: proportion of varation
#' explained by PCA components
#' 
#' @export
#' @examples
#' data(iris)
#' nipals(iris[, 1:4], a = 2)

nipals <- function(X, a, it = 10, tolerance = 1e-4){
    
    if (a > ncol(X)) {
	stop("cannot calculate more components than variables", call. = FALSE)
    }

    Xh <- scale(X, center = TRUE, scale = FALSE) #mean-centering of data matrix X
    nr <- 0
    T <- NULL
    P <- NULL
    pcvar <- NULL
    varTotal <- sum(diag(var(Xh)))
    currVar <- varTotal
    precision <- c()
    for (h in 1:a) {
	th <- Xh[,1] # starting value for th is 1st column of Xh
	ende <- FALSE
	# 3 inner steps of NIPALS algorithm
	while (!ende){
	    nr <- nr + 1
	
	    # the result of matrix multiplication operation (%*%) is a matrix of a single
	    # valule. A matrix cannot multiply another using scalar multiplication (*).
	    # as.vector convert a value of class matrix to a value of class double.
	    # (A'*B)' = B'*A
	    ph <- t((t(th) %*% Xh) * as.vector(1/(t(th) %*% th)))	#LS regression for ph
	    ph <- ph * as.vector(1/sqrt(t(ph) %*% ph))		#normalization of ph
	    thnew <- t(t(ph) %*% t(Xh) * as.vector(1 / (t(ph) %*% ph)))	#LS regression for th
	    prec <- t(th - thnew) %*% (th - thnew)	# calculate precision
	    precision <- append(precision, prec) # record precision
	    th <- thnew	#refresh th in any case
	    #check convergence of th
	    if (prec <= (tolerance ^ 2)) {
		ende <- TRUE
	    } else if (it <= nr) { #too many iterations
		ende <- TRUE
		cat("\nWARNING! Iteration stop in h=",h," without convergence!\n\n")
		}
	    }
	Xh <- Xh - (th %*% t(ph))	#calculate new Xh
	T <- cbind(T,th)	#build matrix T
	P <- cbind(P,ph)	#build matrix P
	oldVar <- currVar
	currVar <- sum(diag(var(Xh)))
	pcvar <- c(pcvar,(oldVar - currVar) / varTotal)
	nr <- 0
    }
    list(T = T, P = P, pcvar = pcvar, precision = precision)
}

