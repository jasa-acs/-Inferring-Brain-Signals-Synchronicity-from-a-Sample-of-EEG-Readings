#' MIC_sim Helper: pars
#'
#' \code{pars} determines the AR(2) coefficient with desired oscillation properties: peak location, peak width and
#'  sampling frequency.
#'
#' @param eta numeric, peak location
#' @param M   numeric, dispersion of power with narrower peak as M approximates 1, default at 1.1
#' @param fs  sampling frequency
#' @return vector of AR(2) coefficients
#' @seealso \code{\link{MIC_sim}} for its usage
#'

pars<- function(eta, M = 1.1, fs){
  phi1 <- - 1 / M ^ 2
  phi2 <- 2 * cos(2 * pi * eta / fs) / M
  return(c(phi2, phi1))
}

#' Estimate Spectral Density using Parzen Lag Window
#'
#' \code{spec.parzen} calculate the spectral estimate based on a Fourier transform of a
#'   truncated and Parzen window smoothed auto-covariance function (ACF).
#'
#' The raw periodogram \code{spec.pgram} is not a consistent estimator of the spectral density,
#'   therefore a class of lag window estimators are considered as surrogates in practice achieving
#'   smootheness and consistency.
#'
#' Parzen window estimator works the best when the true spectrum is continuous, and specifically
#'   have peak concentrations at certain frequencies. Such estimators operates on times series
#'   that are presumably zero-mean and stationary, so that demean and detrending are highly
#'   recommended on \code{x} before implementation.
#'
#' @param x vector, a univariate time series.
#' @param a integer, max lag to truncate the sample ACF estimates, default \code{length(x)-1} (no truncation).
#'   It has to be smaller than the length fo time seires \code{x}.
#' @param nn integer, resolution of the spectral estimates (number of points on frequency domain), default at 512
#'
#' @return List of the following items:
#'   \item{\code{freq}}{frequencies where PSD is evaluated at}
#'   \item{\code{spec}}{estimated correlogram with truncation \code{a} and Parzen window smoothed ACF}
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' spec <- spec.parzen(x, a = 50, nn = 50)
#'
#' plot(spec.pgram(x, plot = "F")$spec, ylab = 'spec')
#' lines(spec$spec)
#' }
#' @export

spec.parzen <- function( x,
                         a  = length(x)-1,
                         nn = 512){
  # Detrrend and Demean:
  x <- lm(x~c(1:length(x)))$residuals
  # tapering
  x <- spec.taper(x, p=0.1)
  # Compute the smoothed periodogram using a Parzen window
  if (a >= length(x)){
    warning("The lag a is too big, reset to 'length(x)-1'.")
    a = length(x) - 1
  }
  freq = seq(0,0.5,length.out = nn+1)[-1]
  spec = specParzen(ts = x, lag = a, maxf = floor(length(x)/2), outn = nn)
  return(list(freq = freq, spec = spec))
}

#' Cluster Aligner
#'
#' \code{align} aligns two cluster assignments by a maximally concordant relabeling, implemented
#'   based on Cpp function \code{\link{clustalign}}
#'
#' This aligner can take both \emph{vector} and \emph{matrix} clustering assignments, and aligns
#'   the second assignment against the first assignment. It reshuffles the second assignment labels
#'   such that two assignemnts are maximally matched. This procedure facilitates a direct and
#'   fair comparision between two clustering result.
#'
#' @param refL vector or matrix, reference label.
#' @param L vector or matrix, label to be aligned.
#' @param type what type of input format: "\code{vec}" or "\code{mat}", default as 'vec'.
#' @return \code{L}, aligned against \code{refL} with maximal concordance.
#'
#' @examples
#' ref <- c(1:5)
#' L <- c(5:1)
#'
#' align(ref, L, type = 'vec')
#'
#' @export

align <- function(refL, L, type = 'vec'){
  # if(type != 'vec' & type != 'mat') stop("Please input type as vec or mat")
  if(!is.numeric(refL) | !is.numeric(L)) stop("refL and L needs to be numeric vector or matrix")
  inType = pmatch(type,c("vec", "mat"), nomatch = 0)
  if (inType == 1){
    maxk = max(refL, L)
    m1 = outer(c(1:maxk), refL, FUN = "==")+0 ; m2 = outer(c(1:maxk), L, FUN = "==")+0
    clustalign(m2, m1)
    L = as.vector(c(1:maxk) %*% m2)
  } else if (inType == 2) {
    clustalign(L,refL)
  } else {
    stop("Input type 'type' must be 'vec' or 'mat'.")
  }
  return(L)
}
