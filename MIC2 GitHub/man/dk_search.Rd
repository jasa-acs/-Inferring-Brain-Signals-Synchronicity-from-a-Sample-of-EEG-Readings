% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.R
\name{dk_search}
\alias{dk_search}
\title{d-K Search Algorithm}
\usage{
dk_search(X_array, max_d = 10, n.iter, par.spec = c(100, 50, 100),
  par.win = c(1, 0), unit_len = T)
}
\arguments{
\item{X_array}{list of data arrays, each of which is organized as No.objects * No.observations * No.segments.}

\item{max_d}{integer, maximal value of \code{d} and \code{K} to be considered}

\item{n.iter}{integer, number of iterations in \code{\link{MIC}} fitting}

\item{par.spec}{vector, spectral estimation parameters the same as \code{par.spec} in \code{\link{MIC_prep}}}

\item{par.win}{vector, epoch smoothing parameters see \code{\link{MIC_prep}} and \code{\link{SpecSim}}}

\item{unit_len}{boolean, whether to use normalized eigen-Laplacian in \code{\link{MIC_prep}}}
}
\value{
A matrix of searching trajectory
  \item{\code{d}}{dimensionality}
  \item{\code{K}}{prespecified number of clusters}
  \item{\code{BIC}}{BIC accordingly}
  \item{\code{COH}}{Adjusted coherence accordingly}
}
\description{
\code{dk_search} selects the optimal eigen-Laplacian dimensionality (\code{d}) and the number of clusters (\code{K})
  jointly based on model assessment criteria like Bayesian Information Criterion (BIC) and adjusted adherence.
}
\details{
Details of the procedure is available in our manuscript.
}
\examples{
\dontrun{
# An example here

## Time series simulation:
  sim <- MIC_sim(alpha = 0.9,  nsub = 10, fs = 200, segs = 10)


## d,K searching: \\strong{(approx. 3 mins running time)}
  dk_search(sim$Data, max_d = 10, n.iter = 10000, par.spec = c(50,50,256), par.win = c(4,2))

}

}
\references{
Qian Li, Damla Senturk, Catherine A. Sugar, Shanali Jeste, Charlotte DiStefano, Joel Frohlich, Donatello Telesca
  "\emph{Inferring Brain Signals Synchronicity from a Sample of EEG Readings}".
}
