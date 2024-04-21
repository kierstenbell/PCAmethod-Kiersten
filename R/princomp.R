#  princomp
#' @title An interactive principle component analysis function
#' @description This function allows you to ...
#'
#' @param x a numeric matrix or data frame with data for PCA
#' @param cor a logical value indicating whether the calculation should use the correlation matrix instead of the covariance matrix. Defaults to FALSE.
#' @param scores a logical value indicating whether teh score on each principal component should be calculated. Defaults to TRUE.
#' @param covmat a covariance matrix or list as returned by `cov.wt`. If TRUE, this is used rather than the covariance matrix of x.
#' @param fix_sign a logical value indicating whether or not the signs of the loadings and scores should be chosen so that the first element of each loading is non-negative.
#' @param formula a formula with no response variable, referring only to numeric variables
#' @param data an optional data frame containing the variables in the formula `formula`. By default the variables are taken from environment(formula)
#' @param subset an optional vector used to select rows (observations) of the data matrix `x`
#' @param na.action a function which indicates what should happen when the data contains NAs. The default is set by the `na.action` setting of `options` and is `na.fail` if that is unset. The 'factory-fresh' default is `na.omit`.
#' @param ... arguments passed to or from other methods. If `x` is a formula one might specify `cor` or `scores`
#' @param object Object of class inheriting from "princomp"
#' @newdata an optional data frame or matrix in which to look for variables with which to predict. If omitted, the scores are used. If the original fit used a formula or a data frame ora a matrix with column names, `newdata` must contain columns with the same names. Otherwise, it must contain the same number of columns to be used in the same order.
# @keywords
# @export
# @examples
#

princomp <- function(x, ...) UseMethod("princomp")

## use formula to allow update() to be used.
princomp.formula <- function(formula, data = NULL, subset, na.action, ...)
{
  mt <- terms(formula, data = data)
  if(attr(mt, "response") > 0) stop("response not allowed in formula")
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  mf$... <- NULL
  ## need stats:: for non-standard evaluation
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval.parent(mf)
  ## this is not a `standard' model-fitting function,
  ## so no need to consider contrasts or levels
  if (.check_vars_numeric(mf))
    stop("PCA applies only to numerical variables")
  na.act <- attr(mf, "na.action")
  mt <- attr(mf, "terms") # allow model.frame to update it
  attr(mt, "intercept") <- 0
  x <- model.matrix(mt, mf)
  res <- princomp.default(x, ...)
  ## fix up call to refer to the generic, but leave arg name as `formula'
  cl[[1L]] <- as.name("princomp")
  res$call <- cl
  if(!is.null(na.act)) {
    res$na.action <- na.act # not currently used
    if(!is.null(sc <- res$scores))
      res$scores <- napredict(na.act, sc)
  }
  res
}

princomp.default <-
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
    ## check if cor parameter selected
    if (cor) {
      sds <- sqrt(diag(cv)) ## standard deviations?
      if(any(sds == 0))
        stop("cannot use 'cor = TRUE' with a constant variable")
      cv <- cv/(sds %o% sds) ## %o% is the outer product of arrays
    }
    edc <- eigen(cv, symmetric = TRUE)
    ev <- edc$values
    if (any(neg <- ev < 0)) { # S-PLUS sets all := 0
      ## 9 * : on Solaris found case where 5.59 was needed (MM)
      if (any(ev[neg] < - 9 * .Machine$double.eps * ev[1L]))
        stop("covariance matrix is not non-negative definite")
      else
        ev[neg] <- 0
    }
    cn <- paste0("Comp.", 1L:ncol(cv))
    names(ev) <- cn
    dimnames(edc$vectors) <- if(missing(x))
      list(dimnames(cv)[[2L]], cn) else list(dimnames(x)[[2L]], cn)
    sdev <- sqrt(ev)
    sc <- setNames(if (cor) sds else rep.int(1, ncol(cv)),
                   colnames(cv))
    scr <- if (scores && !missing(x) && !is.null(cen))
      scale(z, center = cen, scale = sc) %*% edc$vectors
    if (is.null(cen)) cen <- rep(NA_real_, nrow(cv))
    edc <- list(sdev = sdev,
                loadings = structure(edc$vectors, class="loadings"),
                center = cen, scale = sc, n.obs = n.obs,
                scores = scr, call = cl)
    ## The Splus function also return list elements factor.sdev,
    ## correlations and coef, but these are not documented in the help.
    ## coef seems to equal load.  The Splus function also returns list
    ## element terms which is not supported here.
    class(edc) <- "princomp"
    edc
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
