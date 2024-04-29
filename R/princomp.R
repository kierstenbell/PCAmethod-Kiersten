#  princomp
#' @title An interactive principle component analysis function
#' @description This function allows you to ...
#'
#' @param x a numeric matrix or data frame with data for PCA
#' @param cor a logical value indicating whether the calculation should use the correlation matrix instead of the covariance matrix. Can only use if there are no constant variables. Defaults to FALSE.
#' @param scores a logical value indicating whether the score on each principal component should be calculated. Defaults to TRUE.
#' @param covmat a covariance matrix or list as returned by `cov.wt`. If TRUE, this is used rather than the covariance matrix of x.
#' @param fix_sign a logical value indicating whether or not the signs of the loadings and scores should be chosen so that the first element of each loading is non-negative.
#' @param subset an optional vector used to select rows (observations) of the data matrix `x`
#' @keywords pca
#' @export
#' @examples
#' princomp(x = music_clean)

princomp <-
  function(x, cor = FALSE, scores = TRUE, covmat = NULL,
           subset = rep_len(TRUE, nrow(as.matrix(x))), ...) #default subset is the whole matrix.
  {
    chkDots(...)
    cl <- match.call()
    cl[[1L]] <- as.name("princomp")
    ## coerce x into a matrix. Uses a default subset of whole matrix, unless otherwise specified in input.
    z <- if(!missing(x)) as.matrix(x)[subset, , drop = FALSE]
    ## check if covmat exists as a list.
    if (is.list(covmat)) {
      if(any(is.na(match(c("cov", "n.obs"), names(covmat)))))
        stop("'covmat' is not a valid covariance list")
      cv <- covmat$cov
      n.obs <- covmat$n.obs
      cen <- covmat$center
    ## check if covmat exists as a matrix
    } else if(is.matrix(covmat)) {
      if(!missing(x)) ## warn only here; x is used for scores when we have 'cen'
        warning("both 'x' and 'covmat' were supplied: 'x' will be ignored")
      cv <- covmat
      n.obs <- NA
      cen <- NULL
    ## continue if covmat variable is NULL (the default)
    } else if(is.null(covmat)){
      dn <- dim(z) ## get dimensions of matrix z
      ## check if matrix is capable of running a PCA
      if(dn[1L] < dn[2L])
        stop("'princomp' can only be used with more units than variables")
      ## get covariance matrix from input z
      covmat <- cov.wt(z)             # returns list, cov() does not
      ## store the number of observations
      n.obs <- covmat$n.obs
      ## normalize covariance data
      cv <- covmat$cov * (1 - 1/n.obs)# for S-PLUS compatibility
      cen <- covmat$center ## get centers
    } else stop("'covmat' is of unknown type") ## break if covmat input is neither null, list, nor matrix
    if(!is.numeric(cv)) stop("PCA applies only to numerical variables")
    ## check if cor parameter selected.
    if (cor) {
      sds <- sqrt(diag(cv)) ## standard deviations of covariance??
      ## cancels operation if there is a constant variable (sds=0) in dataset. Correlation cannot be used here.
      if(any(sds == 0))
        stop("cannot use 'cor = TRUE' with a constant variable")
      cv <- cv/(sds %o% sds) ## calculate covariance matrix. %o% is the outer product of arrays
    }
    ## computes the eigenvalues and eigenvectors of numeric or complex matrices
    ## inputs the covariance matrix `cv` calculated above
    edc <- eigen(cv, symmetric = TRUE)
    ## extracts the eigenvalues
    ev <- edc$values
    if (any(neg <- ev < 0)) { # S-PLUS sets all := 0
      ## 9 * : on Solaris found case where 5.59 was needed (MM)
      if (any(ev[neg] < - 9 * .Machine$double.eps * ev[1L]))
        stop("covariance matrix is not non-negative definite")
      else
        ev[neg] <- 0
    }
    ## store the total number of principal components based on number of columns in covariance matrix
    cn <- paste0("Comp.", 1L:ncol(cv))
    ## renames the eigenvalues with their corresponding principal component number
    names(ev) <- cn
    ## rename eigenvectors if missing x ??
    dimnames(edc$vectors) <- if(missing(x))
      list(dimnames(cv)[[2L]], cn) else list(dimnames(x)[[2L]], cn)
    ## standard deviations of the principal componenets (square root of eigenvalues)
    sdev <- sqrt(ev)
    ## set scaling factor for the
    sc <- setNames(if (cor) sds else rep.int(1, ncol(cv)),
                   colnames(cv))
    ## scale z by the scores if all conditions are TRUE
    ## only will work if scores are given in original input
    scr <- if (scores && !missing(x) && !is.null(cen))
      scale(z, center = cen, scale = sc) %*% edc$vectors  ## %*% multiplies the scaling matrix with the eigenvector matrix
    ## add NAs to null centers
    if (is.null(cen)) cen <- rep(NA_real_, nrow(cv))
    ## recompile the eigenvector (PCA) information
    edc <- list(sdev = sdev,
                loadings = structure(edc$vectors, class="loadings"), ## loadings are the correlations between the component and the original variables; i.e., how much of the variation in a variable is explained by the component
                center = cen, scale = sc, n.obs = n.obs,
                scores = scr, call = cl)
    ## The Splus function also return list elements factor.sdev,
    ## correlations and coef, but these are not documented in the help.
    ## coef seems to equal load.  The Splus function also returns list
    ## element terms which is not supported here.
    class(edc) <- "princomp"
    edc ##final return is the eigenvectors, now renamed
  }

print.princomp <- function(x, ...)
{
  cat("Call:\n"); dput(x$call, control=NULL)
  cat("\nStandard deviations:\n")
  print(x$sdev, ...)
  cat("\n", length(x$scale), " variables and ", x$n.obs,
      "observations.\n")
  invisible(x)
}
