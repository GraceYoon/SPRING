#' SPRING main function
#' Semi-Parametric Rank-based approach for INference in Graphical model. SPRING follows the neighborhood selection methodology outlined in "mb" Meinshausen and Buhlmann and
#'
#' @param data n by p matrix of microbiome data. Each row represents each subject/sample and each column represents each OTU (operational taxonomic unit).
#' @param quantitative default is FALSE, which means input "data" is compositional data, which will be normalized using mclr transformation within a function. If TRUE, it means "quantitative" counts are input and no normalization will be applied.
#' @param method graph estimation methods. Currently, only "mb" method is available.
#' @param lambda.min.ratio default is 0.01
#' @param nlambda default is 50.
#' @param lambdaseq a sequence of decreasing positive numbers to control the regularization. Users can specify a sequence to override the default sequence.
#' @param seed the seed for subsampling
#' @param ncores number of cores to use subsampling
#' @param thresh threshold for StARS selection criterion. 0.1 is recommended.
#' @param subsample.ratio 0.8 is default. The recommended values are 10*sqrt(n)/n for n > 144 or 0.8 otherwise.
#' @param rep.num the repetition number of subsampling
#'
#' @return \code{SPRING} returns a data.frame containing
#' @importFrom huge huge.mb
#' @importFrom pulsar pulsar
#' @importFrom mixedCCA estimateR
#' @import SpiecEasi
#'
#' @export
#'
#' @examples
#' library(SpiecEasi)
#' data(amgut1.filt)
#' fit <- SPRING(data = amgut1.filt[, 1:10], quantitative = FALSE, nlambda = 20, rep.num = 10)
SPRING <- function(data, quantitative = FALSE, method = "mb", lambda.min.ratio = 1e-2, nlambda = 50, lambdaseq = NULL, seed = 10010, ncores = 2, thresh = 0.1, subsample.ratio = 0.8, rep.num = 50){
  p <- ncol(data)
  if (quantitative){
    qdat <- data
  } else {
    qdat <- mclr(data)
  }
  Kcor <- mixedCCA::estimateR(qdat, type = "trunc")$R

  if(is.null(lambdaseq)){
    # generate lambda sequence
    lambda.max <- max(max(Kcor-diag(p)), -min(Kcor-diag(p)))
    lambda.min <- lambda.min.ratio * lambda.max
    lambdaseq <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  }

  if(method == "mb"){
    fun <- hugeKmb
  }

  out1.K_count <- pulsar::pulsar(qdat, fun = fun, fargs = list(lambda = lambdaseq), rep.num = rep.num, criterion = 'stars', seed = seed, ncores = ncores, thresh = thresh, subsample.ratio = subsample.ratio)

  fit1.K_count <- pulsar::refit(out1.K_count)

  return(list(Kcor = Kcor, output = out1.K_count, fit = fit1.K_count, lambdaseq = lambdaseq))
}
