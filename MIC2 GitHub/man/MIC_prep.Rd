% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MIC_prep.R
\name{MIC_prep}
\alias{MIC_prep}
\title{Data preparation for MIC}
\usage{
MIC_prep(X, d, exclude = NULL, par.spec = NULL, par.win = c(1, 0),
  unit_len = FALSE, spec = FALSE)
}
\arguments{
\item{X}{3-d array, organized as No.objects * No.observations * No.segments, see simulated example by \code{\link{MIC_sim}}}

\item{d}{integer, dimensionality of the eigen-Laplacian representation}

\item{exclude}{vector of integers, indicating objects to be excluded}

\item{par.spec}{vector of parameters as \code{spec.lag}, \code{spec.mfreq} and \code{spec.len} as in \code{\link{SpecSim}}}

\item{par.win}{vector, moving average smoothing parameters, see \code{\link{SpecSim}}, default at \code{c(1,0)}}

\item{unit_len}{boolean, to use normalized eigen vector in \code{\link{EigLap}} (default \code{FALSE})}

\item{spec}{boolean, to return spectral estimates or MIC ready data format. (default \code{FALSE})}
}
\value{
List of data matrices, each with No.objects rows and \code{d} columns.
}
\description{
\code{MIC_prep} prepares stuctured data for \code{MIC}. It sequentially carries out signal detrending (linear),
  tapering, epoch smoothing, distance characterization and dimensional reduction.
}
\examples{
\dontrun{
# Simulated data:
ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 100)

# Data preparation on subject 1
sub1 <- MIC_prep(ts_sim$Data[[1]], d = 3, par.spec = c(80,50), par.win = c(3, 1))


# Data structure: D / No.channels / No.Epochs
dim(sub1)


# To visualize preprocessed data on Epoch 1
ep1 <- t(sub1[,,1])
plot(ep1[,c(1,2)], col = ts_sim$Ci[,1], xlab = 'dim1', ylab = 'dim2')
}
}
\seealso{
\code{\link{MIC}} for its usage and \code{\link{MIC_sim}} for time series simulation.
}
