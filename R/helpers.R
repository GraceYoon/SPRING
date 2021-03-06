# Rfunctions for simple example for SPRING method
# Yoon, Gaynanova and M\"{u}eller (2019) Frontiers in Genetics, Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data.
# doi:10.3389/fgene.2019.00516

# To implement Kendall correlation estimates on huge.
# original huge function can only take covariance or data matrix.
# huge.mb funciton returns beta values (MB coefficient estimates)
# huge function with method="mb" does not return beta values.
# to do network visualization (for edge color, I need beta)


#' Internal wrapper function to implement rank-based correlation to huge.mb function in "huge" package.
#'
#' @param data n by p matrix data. usually through pulsar, data will receive subsamples.
#' @param lambda a vector of lambda values
#' @param type a type of variables. "trunc" is default.
#' @param sym "or" is the symmetrizing rule of the output graphs. If sym = "and", the edge between node i and node j is selected ONLY when both node i and node j are selected as neighbors for each other. If sym = "or", the edge is selected when either node i or node j is selected as the neighbor for each other. The default value is "or". (refer to huge manual)
#' @param verbose If \code{verbose = FALSE}, tracing information printing for HUGE (High-dimensional Undirected Graph Estimation) with a specified method (currently "mb" is only available) is disabled. The default value is TRUE.
#' @param verboseR If \code{verboseR = FALSE}, printing information whetehr nearPD is used or not is disabled. The defalut value is TRUE.
#' @param Rmethod The calculation method of latent correlation. Either "original" method or "approx". If \code{Rmethod = "approx"}, multilinear approximation method is used, which is much faster than the original method. If \code{Rmethod = "original"}, optimization of the bridge inverse function is used. The default is "approx".
#' @param tol Desired accuracy when calculating the solution of bridge function in estimateR function.
#'
#' @return \code{hugeKmb} returns a data.frame containing
#' \itemize{
#'      \item{beta: }
#'      \item{path: }{a list of}
#'      \item{df: }
#' }
#'
#' @importFrom huge huge.mb
#' @export
#'
hugeKmb <- function(data, lambda, type = "trunc", sym = "or", verbose = TRUE, verboseR = TRUE, Rmethod = "approx", tol = 1e-6) {
  S    <- mixedCCA::estimateR(data, type = type, method = Rmethod, tol = tol, verbose = verboseR)$R
  est  <- huge::huge.mb(S, lambda, sym = sym, verbose = verbose)
  est
}




#' Modified central log ratio (mclr) transformation
#'
#' @param dat raw count data or compositional data (n by p) does not matter.
#' @param base exp(1) for natural log
#' @param tol tolerance for checking zeros

# For eps and atleast, users do not have to specify any values. Default should be enough.
#' @param eps epsilon in eq (2) of the paper "Yoon, Gaynanova, M\"{u}ller (2019), Frontiers in Genetics". positive shifts to all non-zero compositions. Refer to the paper for more details. eps = absolute value of minimum of log ratio counts plus c.
#' @param atleast default value is 1. Constant c which ensures all nonzero values to be strictly positive. default is 1.
#'
#'
#' @return \code{mclr} returns a data matrix of the same dimension with input data matrix.
#' @export
#'
#' @examples
#' data(QMP)
#' RMP <- QMP/rowSums(QMP)
#' mclr_RMP <- mclr(RMP)
#'
mclr <- function(dat, base = exp(1), tol = 1e-16, eps = NULL, atleast = 1){
  dat <- as.matrix(dat)
  nzero <- (dat >= tol)  # index for nonzero part
  LOG <- ifelse(nzero, log(dat, base), 0.0) # take log for only nonzero values. zeros stay as zeros.

  # centralize by the log of "geometric mean of only nonzero part" # it should be calculated by each row.
  if (nrow(dat) > 1){
    clrdat <- ifelse(nzero, LOG - rowMeans(LOG)/rowMeans(nzero), 0.0)
  } else if (nrow(dat) == 1){
    clrdat <- ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
  }

  if (is.null(eps)){
    if(atleast < 0){
      warning("atleast should be positive. The functions uses default value 1 instead.")
      atleast = 1
    }
    if( min(clrdat) < 0 ){ # to find the smallest negative value and add 1 to shift all data larger than zero.
      positivecst <- abs(min(clrdat)) + atleast # "atleast" has default 1.
    }else{
      positivecst <- 0
    }
    # positive shift
    ADDpos <- ifelse(nzero, clrdat + positivecst, 0.0) ## make all non-zero values strictly positive.
    return(ADDpos)
  } else if(eps == 0){
    ## no shift. clr transform applied to non-zero proportions only. without pseudo count.
    return(clrdat)
  } else if(eps > 0){
    ## use user-defined eps for additional positive shift.
    ADDpos <- ifelse(nzero, clrdat + eps, 0.0)
    return(ADDpos)
  } else {
    stop("check your eps value for additional positive shift. Otherwise, leave it as NULL.")
  }
}
