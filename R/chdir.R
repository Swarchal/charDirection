#' Characteristic Direction
#'
#' Calculates the characteristic direction for a phenotypic dataset
#'
#' @param ctrl matrix of control data
#' @param expm matrix of experiment (perturbed) data
#' @param samples vector of drug names
#' @param r numeric between 0 and 1. Regularised term, a parameter that smooths
#'	the covariance matrix and reduces potential noise in the dataset. The
#'	default value is 0, no regularisation.
#'
#' @return vector of n-components representing the characteristic direction of
#'	the dataset. \code{n} equals the number of features in the dataset.
#'	b is also a matrix object, and is sorted by its components' absolute
#'	values in descending order.
#' @export

chdir <- function(ctrl, expm, samples, r = 1){
    
    if (dim(ctrl)[1] != dim(expm)[1]){
    	stop('Control expression data must have equal number of genes as experiment expression data!')
    }
    
    if (any(is.na(ctrl)) || any(is.na(expm))){
    	stop('Control expression data and experiment expression data have to be real numbers. NA was found!')
    }
    
    
    # There should be variance in expression values of each gene. If  
    # gene expression values of a gene are constant, it would dramatically
    # affect the LDA caculation and results in a wrong answer.
    constantThreshold <- 1e-5
    ctrlConstantGenes <- diag(var(t(ctrl))) < constantThreshold
    expmConstantGenes <- diag(var(t(expm))) < constantThreshold
    
    if (any(ctrlConstantGenes)){
    	errMes <- sprintf('%s row(s) in control expression data are constant. Consider Removing the row(s).',
		paste(as.character(which(ctrlConstantGenes)),collapse=','))
    	stop(errMes, call. = FALSE)
    }else if(any(expmConstantGenes)){
    	errMes <- sprintf('%s row(s) in experiment expression data are constant. Consider Removing the row(s).',
		paste(as.character(which(expmConstantGenes)),collapse=','))
    	stop(errMes, call. = FALSE)
    }
    
    # place control gene expression data and experiment gene expression data into
    # one matrix
    combinedData <- cbind(ctrl, expm)
    
    # get the number of samples, namely, the total number of replicates in  control 
    # and experiment. 
    samplesCount <- dim(combinedData)[2]
    
    # the number of output components desired from PCA. We only want to calculate
    # the chdir in a subspace that capture most variance in order to save computation 
    # workload. The number is set 20 because considering the number of genes usually 
    # present in an expression matrix 20 components would capture most of the variance.
    componentsCount <- min(c(samplesCount - 1, 20))
    
    
    # use the nipals PCA algorithm to calculate R, V, and pcvars. nipals algorithm
    # has better performance than the algorithm used by R's builtin PCA function.
    # R are scores and V are coefficients or loadings. pcvars are the variances 
    # captured by each component 
    pcaRes <- nipals(t(combinedData), componentsCount, 1e5)
    R <- pcaRes$T # PCA scores
    V <- pcaRes$P # PCA loadings
    pcvars <- pcaRes$pcvar
    
    
    # we only want components that cpature 95% of the total variance or a little above.
    # cutIdx is the index of the compoenent, within which the variance is just equal
    # to or a little greater than 95% of the total.
    cutIdx <- which(cumsum(pcvars) > 0.95)
    if (length(cutIdx)==0){
    	cutIdx <- componentsCount
    } else{
    	cutIdx <- cutIdx[1]
    }
    
    # slice R and V to only that number of components.
    R <- R[,1:cutIdx]
    V <- V[,1:cutIdx]
    
    # the difference between experiment mean and control mean.
    meanvec <- rowMeans(expm) - rowMeans(ctrl) 
    
    # shruken covariance matrix as denominator for LDA
    shrunkMats <- shrink_mat(R, samplesCount, r)

    # The LDA formula.
    #  V%*%solve(shrunkMats)%*%t(V) transforms the covariance matrix from the subspace to full space.
    b <- V %*% solve(shrunkMats) %*% t(V) %*% meanvec
    
    # normlize b to unit vector
    b <- b * as.vector(sqrt(1 / t(b) %*% b))
    
    # sort b to by its components' absolute value in decreasing order and get the 
    # sort index
    sortRes <- sort(abs(b), decreasing = TRUE, index.return = TRUE)
    
    # sort b by the sort index
    bSorted <- as.matrix(b[sortRes$ix])
    # sort genes by the sort index
    samplesSorted <- samples[sortRes$ix]
    # assign genesSorted as the row names of bSorted
    rownames(bSorted) <- samplesSorted
    return(bSorted)
}

