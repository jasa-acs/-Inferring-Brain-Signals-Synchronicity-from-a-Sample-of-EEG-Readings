#' Data preparation for MIC
#'
#' \code{MIC_prep} prepares stuctured data for \code{MIC}. It sequentially carries out signal detrending (linear),
#'   tapering, epoch smoothing, distance characterization and dimensional reduction.
#'
#' @param X 3-d array, organized as No.objects * No.observations * No.segments, see simulated example by \code{\link{MIC_sim}}
#' @param d integer, dimensionality of the eigen-Laplacian representation
#' @param exclude vector of integers, indicating objects to be excluded
#' @param par.spec vector of parameters as \code{spec.lag}, \code{spec.mfreq} and \code{spec.len} as in \code{\link{SpecSim}}
#' @param par.win vector, moving average smoothing parameters, see \code{\link{SpecSim}}, default at \code{c(1,0)}
#' @param unit_len boolean, to use normalized eigen vector in \code{\link{EigLap}} (default \code{FALSE})
#' @param spec boolean, to return spectral estimates or MIC ready data format. (default \code{FALSE})
#' @return List of data matrices, each with No.objects rows and \code{d} columns.
#' @examples
#' \dontrun{
#' # Simulated data:
#' ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 100)
#'
#' # Data preparation on subject 1
#' sub1 <- MIC_prep(ts_sim$Data[[1]], d = 3, par.spec = c(80,50), par.win = c(3, 1))
#'
#'
#' # Data structure: D / No.channels / No.Epochs
#' dim(sub1)
#'
#'
#' # To visualize preprocessed data on Epoch 1
#' ep1 <- t(sub1[,,1])
#' plot(ep1[,c(1,2)], col = ts_sim$Ci[,1], xlab = 'dim1', ylab = 'dim2')
#' }
#'@seealso \code{\link{MIC}} for its usage and \code{\link{MIC_sim}} for time series simulation.
#'
#' @export
MIC_prep <- function(X, d,
                     exclude = NULL,
                     par.spec= NULL,
                     par.win = c(1, 0),
                     unit_len  = FALSE,
                     spec  = FALSE){
  if (!is.null(exclude)) X <- X[- exclude, , ]         #chanel exclusion
  segs <- dim(X) [3]; fs <- dim(X) [2]; nc <- dim(X) [1]  #extract consts
  # Detrend (Linear) and tapering with .10 on both end
  matx <- cbind(1, c(1:fs)/fs)
  X    <- apply(X, c(1,3), function(v) spec.taper(x = (v - matx %*% qr.solve(matx,v)), p = 0.1))
  # Spectral estimation parameters (default)
  spec.lag = fs-1; spec.mfreq = floor(fs/2); spec.len = 512
  if(length(par.spec) > 3) {
    stop("'par.spec' has 3 parameters at maximum.")
  } else if (length(par.spec)==3) {
    spec.lag=par.spec[1];spec.mfreq=par.spec[2];spec.len=par.spec[3];
  } else if (length(par.spec)==2) {
    spec.lag=par.spec[1];spec.mfreq=par.spec[2];
  } else if (length(par.spec)==1) {
    spec.lag=par.spec[1]
  }
  if(!spec){
    # Cpp based spectral similarity
    ep_sim  <- SpecSim(ts = X, lag = spec.lag, wn = spec.mfreq, win = par.win[1], overlap = par.win[2], specN = spec.len)
    # Eigen Laplacian projection
    ep_eig  <- EigLap(data = ep_sim, D = d, normal = unit_len)
  } else {
    ep_eig <- SpecOnly(ts = X, lag = spec.lag, wn = spec.mfreq, win = par.win[1], overlap = par.win[2], specN = spec.len)
  }
  return(ep_eig)
}
