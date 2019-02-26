#' 128-EEGNet channel plots
#'
#' \code{EEGplot} plots the clustering results of 124 channels on a 2D projected scalp plot.
#'
#' @param clust vector, cluster labels for 124 channels.
#' @param color vector, colors for K clusters in the plot.
#' @examples
#' \dontrun{
#' # The channel plot:
#'
#' EEGplot(rep(1,124))
#'
#' }
#'
#' @importFrom graphics box plot points
#' @export
EEGplot <- function(clust, color = NULL)
{ # A few checks
  K <- max(clust); if(length(clust)!=124) stop("Length of clust needs to be 124!")
  if(is.null(color)) color <- rep('black', K)
  if(length(color)<K) stop("Number of colors needs to be greater than K (no.cluster)!")
  # Initiate
  plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='', xlim = c(-1,1),ylim=c(-1,1))
  # Box
  box()
  # Head and outer circle
  angs <- seq(0,2*pi,length.out = 500);
  outr <- 1.08; inr <- 0.7; nose <- c(inr,inr+.12,inr)
  points(x=outr*sin(angs),y=outr*cos(angs), type='l', col = 'grey')
  points(x=inr*sin(angs),y=inr*cos(angs), type='l', col = 'black', lwd=3)
  # Nose
  noseloc <- c(-0.12,0,0.12)
  points(x=noseloc,y=sqrt(nose^2-noseloc^2), type='l', col='black', lwd=3)
  if(length(color) > 1){
    # Region Coloring
    samp_points <- expand.grid(x = seq(-1.08,1.08,by=.01),
                               y = seq(-1.08,1.08,by=.01))
    samp_lab <- (samp_points$x^2+samp_points$y^2) < 1.08^2
    samp_points <- samp_points[samp_lab,]
    samp_col <- class::knn(EEGcoords[,c(2,3)], samp_points,clust)
    points(x = samp_points$x, y = samp_points$y, pch = ".", col = color[samp_col])
  }
  # channels
  points(x = EEGcoords$x, y = EEGcoords$y, pch = 19, col = color[clust])
}
