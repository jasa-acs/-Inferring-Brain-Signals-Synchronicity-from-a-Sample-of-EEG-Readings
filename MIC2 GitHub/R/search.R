#' d-K Search Algorithm
#'
#' \code{dk_search} selects the optimal eigen-Laplacian dimensionality (\code{d}) and the number of clusters (\code{K})
#'   jointly based on model assessment criteria like Bayesian Information Criterion (BIC) and adjusted adherence.
#'
#' Details of the procedure is available in our manuscript.
#'
#' @references Qian Li, Damla Senturk, Catherine A. Sugar, Shanali Jeste, Charlotte DiStefano, Joel Frohlich, Donatello Telesca
#'   "\emph{Inferring Brain Signals Synchronicity from a Sample of EEG Readings}".
#' @param X_array list of data arrays, each of which is organized as No.objects * No.observations * No.segments.
#' @param max_d integer, maximal value of \code{d} and \code{K} to be considered
#' @param n.iter integer, number of iterations in \code{\link{MIC}} fitting
#' @param par.spec vector, spectral estimation parameters the same as \code{par.spec} in \code{\link{MIC_prep}}
#' @param par.win vector, epoch smoothing parameters see \code{\link{MIC_prep}} and \code{\link{SpecSim}}
#' @param unit_len boolean, whether to use normalized eigen-Laplacian in \code{\link{MIC_prep}}
#' @return A matrix of searching trajectory
#'   \item{\code{d}}{dimensionality}
#'   \item{\code{K}}{prespecified number of clusters}
#'   \item{\code{BIC}}{BIC accordingly}
#'   \item{\code{COH}}{Adjusted coherence accordingly}
#' @examples
#' \dontrun{
#' # An example here
#'
#' ## Time series simulation:
#'   sim <- MIC_sim(alpha = 0.9,  nsub = 10, fs = 200, segs = 10)
#'
#'
#' ## d,K searching: \strong{(approx. 3 mins running time)}
#'   dk_search(sim$Data, max_d = 10, n.iter = 10000, par.spec = c(50,50,256), par.win = c(4,2))
#'
#' }
#'
#' @export
dk_search <- function(X_array,
                      max_d = 10,
                      n.iter,
                      par.spec = c(100,50,100),
                      par.win  = c(1, 0),
                      unit_len = T){
  result <- c()
  # Initiate
  D <- 2
  K <- 2
  while (D <= max_d){
    list_data <- lapply(X_array, function(x) MIC_prep(x, d=D,
                                      par.spec = par.spec, par.win = par.win, unit_len = unit_len))
    if (D > 2) K <- K - 1 # conservative step-back
    output <- MIC(data = list_data, K = K, nit = n.iter, drop = T)
    current_ICs <- as.vector(output$ICs)
    current_BIC <- output$ICs[[1]]
    current_COH <- output$ICs[[2]]
    result <- rbind(result, c(D, K, current_ICs)) # track the trajectory
    rm(output)

    while(K < max_d){
      output <- MIC(list_data, K = K + 1, nit = n.iter)
      result <- rbind(result, c(D, K + 1, as.vector(output$ICs)))
      if (output$ICs[[1]]> current_BIC){
        current_BIC <- output$ICs[[1]]
        current_COH <- output$ICs[[2]]
        K <- K + 1
      } else {
        cat('D = ',D,': best_K = ',K,' COH = ',current_COH,'\n', sep = '')
        cat('-----------------------------------------','\n')
        break
      }
    }
    if(K == max_d){
      warning('Max D has reached!'); break
    }
    if (D == K) break; # converge
    D <- K # next D
    if (D %in% result[,1]) break;
  }
  # if(D != K){
  #   # Then chose D,K by max coherence, locally
  #   D_cand <- result[, 1] >= (D - 1)
  #   cand_result <- result[D_cand, ]
  #   pick <- which.max(cand_result[, 4])
  #   D <- cand_result[pick,1]
  #   K <- cand_result[pick,2]
  # }
  rm(output)
  colnames(result) <- c('d','K','BIC','COH')
  return(result)
}
