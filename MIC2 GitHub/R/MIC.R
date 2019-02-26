#' Multilevel Integrative Clustering
#'
#' \code{MIC} implements integrative clustering on highly structured data. The admissible structure contains repeated measurements
#'   of the same set of units collected from multiple individuals. Clustering is performed integratively therefore interpretable
#'   on all three levels (group/individual/repetition).
#'
#' @param data list of 3d arrays, each of which has the size of d(measurement dims)*p(number of units)*r(number of repetitions).
#'  MIC requires \code{p} the same across arrays. For an admissible format, see \code{\link{MIC_prep}}.
#' @param K integer, number of clusters.
#' @param nit number of iterations for MCMC sampling, with 1/5 used for burn-in.
#' @param thin thining, default at 1 (no thinning).
#' @param drop logical, to remove posterior samples (default TRUE).
#'
#' @return A list of objects with the following components:
#'   \item{\code{clustS}}{group level clustering result}
#'   \item{\code{clustC}}{individual level clustering result}
#'   \item{\code{clustL}}{repetition level clustering result}
#'   \item{\code{mixpr}}{cluster proportions with its posterior summaries}
#'   \item{\code{alpha}}{S-C level adherence and its posterior summaries}
#'   \item{\code{beta}}{C-L level adherence and its posterior summaries}
#'   \item{\code{ICs}}{model assessment measures: BIC and adjusted adherence}
#'
#' @references Qian Li, Damla Senturk, Catherine A. Sugar, Shanali Jeste, Charlotte DiStefano, Joel Frohlich, Donatello Telesca
#'   "\emph{Inferring Brain Signals Synchronicity from a Sample of EEG Readings}".
#' @seealso \code{\link{MIC_prep}} for data preparation, \code{\link{MIC_sim}} for MIC on simulated time series.
#' @examples
#'
#' \dontrun{
#'   # Time Series simulation:
#'   ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 100)
#'
#'
#'   # Data preparation:
#'   list_data <- lapply(ts_sim$Data, function(x)
#'        MIC_prep(X = x, d = 4, par.spec = c(50, 50, 100), par.win = c(3, 1)))
#'
#'
#'   # MIC: (Running time: 11s, 30X faster than version 1)
#'   output <- MIC(data = list_data, K = 4, nit = 10000)
#'
#'
#'
#'   # Clustering accuracy: group level (to be updated)
#'   sum(align(ts_sim$C, output$clustS, type = 'vec') == ts_sim$C) / 40
#'   ## [1] 1
#'
#'
#'   # Clustering accuracy: individual level (to be updated)
#'   true_c <- lapply(1:3,function(i) ts_sim$Ci[,i]);
#'   est_c <- lapply(1:3,function(i) output$clustC[i,]);
#'   mapply(function(x,y) sum(align(x,y,'vec')==x)/40, true_c, est_c);
#'   ## [1] 1 1 1
#' }
#' @export

MIC <- function(data, K, nit, thin = 1, drop=T){
  # Input checks:
  # - Data structure
  ns <- length(data);     np <- dim(data[[1]])[2]
  ne <- c()
  for(i in 1:ns){
    if(length(dim(data[[i]]))!=3)   stop("Data needs to be a list of 3D arrays")
    if(dim(data[[i]])[2]!=np)  stop("Subject units differ, cannot perform MIC")
    ne <- c(ne, dim(data[[i]])[3])
  }
  # - K
  if(round(K) != K)        stop("K needs to be an integer")
  if(K==1)                  stop("No clustering performed when K=1")
  # - Missing values
  ### Advanced settings:
  # - Flexible prior changes
  # - Indiv alpha; beta options
  # - Options to inculde a timer/Display
  #
  # Create Output directory ---------------------------------------------------
  files <- list.files();
  dir   <- files == "post";
  if(length(files[dir]) == 0) dir.create("post");
  #
  ##############################################################################
  ##                      Posterior sample summaries                          ##
  ##############################################################################
  output <- MIC_mcmc(data, K, nit, thin)
  ## Organizing outputs
  ##  - pop cluster
  clustS <- as.vector(c(1:K) %*% output$clustS)
  ##  - individual cluster
  clustC <- t(apply(output$clustC, c(3), function(x) as.vector(c(1:K) %*% x)))
  rownames(clustC) <- paste0('sub',c(1:ns))
  ##  - Read in Ci's
  tab_Ci <- read.table("post/CIs.txt", header = F, sep = " ")
  ##  - Mixing Probability
  mixpr <- c(); mixpr <- cbind(mixpr, output$pi, t(tab_Ci[,1:K]))
  rownames(mixpr) <- paste0("c",c(1:K)); colnames(mixpr) <- c('mean','lo.ci',
                                                              'median','up.ci')
  tab_Ci <- tab_Ci[, -c(1:K)]
  ##  - Alpha's
  alpha <- c(); alpha <- cbind(alpha, output$alpha, t(tab_Ci[,1:ns]))
  rownames(alpha) <- paste0("sub",c(1:ns)); colnames(alpha) <- c('mean','lo.ci',
                                                              'median','up.ci')
  tab_Ci <- tab_Ci[, -c(1:ns)]
  ##  - Beta's
  beta <- list()
  for(sub in 1:ns){
    betam <- c(); betam <- cbind(betam, output$beta[[sub]], t(tab_Ci[,1:ne[sub]]))
    rownames(betam)<- paste0("ep",c(1:ne[sub]))
    colnames(betam) <- c('mean','lo.ci','median','up.ci')
    tab_Ci <- tab_Ci[, -c(1:ne[sub])]; beta[[sub]] <- betam
  }
  ##  - epoch clusters
  clustL <- list()
  for(sub in 1:ns){
    clust <- t(apply(output$clustL[[sub]], c(3), function(x) as.vector(c(1:K) %*% x)))
    rownames(clust) <- paste0('ep',c(1:ne[sub]))
    clustL[[sub]] <- clust
  }
  ##  - IC's
  ICs <- as.vector(output$ICs)
  names(ICs) <- c("BIC", "COH")
  ###################CLEAR MEMORY and DELETE POST###############################
  rm(output, tab_Ci, clust, betam)
  if(drop) unlink("post", recursive = T)
  ## Return
  return(list(
    clustS = clustS,
    clustC = clustC,
    clustL = clustL,
    mixpr  = mixpr,
    alpha  = alpha,
    beta   = beta,
    ICs    = ICs
  ))
}

#' @useDynLib MIC2
#' @importFrom Rcpp sourceCpp
NULL
