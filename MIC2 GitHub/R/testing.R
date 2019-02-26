# Importation of ourside functions
#' @importFrom stats runif rnorm lm arima.sim var acf rexp spec.taper
#' @importFrom utils tail read.table
NULL
if(F){
  ################## dk search ##################
  set.seed(123)
  sim <- MIC_sim(alpha = 0.9,  nsub = 10, fs = 200, segs = 10)
  start <- Sys.time()
  out <- dk_search(sim$Data, max_d = 10, n.iter = 10000, par.spec = c(100,50,256), par.win=c(4,2))
  out
  Sys.time() - start
  ################### simple fit ################
  for(i in 1:300){
    # set.seed(555)
    ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 200)
    list_data <- lapply(ts_sim$Data, function(x) MIC_prep(X = x, d = 4,
                            par.spec = c(100,50,100), par.win = c(2, 0), unit_len = F))
    # ## plots
    # ep <- list_data[[2]][,,3]; par(mfrow = (c(2,2)));
    # plot((ep[1,]), ep[2,], col = ts_sim$Ci[,2]);plot((ep[1,]), ep[3,], col = ts_sim$Ci[,2])
    # plot((ep[1,]), ep[4,], col = ts_sim$Ci[,2]);plot(ep[2,], ep[3,], col = ts_sim$Ci[,2])

    start <- Sys.time()
    output <- MIC(data = list_data, K = 4, nit = 10000, drop = F); Sys.time() - start

    ## truth and estimated:
    ts_sim$C;align(refL = ts_sim$C, L = output$clustS, type = 'vec')
    t(ts_sim$Ci);for(i in 1:dim(ts_sim$Ci)[2]){
      output$clustC[i,] = align(ts_sim$Ci[,i],output$clustC[i,],type='vec')
    }
    output$clustC
  }
  ################### TS and specs ##############
  ts <- arima.sim(n = 100, model = list(ar = pars(20, fs = 100)))
  ts <- ts - mean(ts)
  spec1 <- Mod(fft(ts))^2/length(ts) ### raw spec
  spec2 <- spec.pgram(ts, plot = F, demean = F, detrend = F)$spec  ### spec.pgram
  spec3 <- spec.parzen(ts, a = 99, nn=50) ### spec.parzen
  matplot(cbind(spec1[2:51], spec2, spec3), type = 'l')
  ################### SpecSim, EigLap and EigLapSph #####
  sim <- MIC_sim(alpha = 1, nsub = 1, Ct = c(1,1,1,2,2,2,3,3,3), fs = 200, segs = 10)$Data[[1]]
  sim <- apply(sim, c(1,3), function(x) lm(x~c(1:200))$residuals)
  test <- SpecSim(sim, 60, 50, 10,0,100)
  for(i in 1:2){
    EigLap(test, 2, (i==2)) -> tplot
    plot(t(tplot[,,1]), col = rep(c(1,2,3),c(3,3,3))) #Consider Cartesan-spherical transform
  }

  specs <- apply(sim, c(2,3), function(x) specParzen(x,100,50,100)/sum(specParzen(x,100,50,100)))
  specs <- apply(specs, c(1,2), mean)
  matplot(specs, type = 'l', col = c(1,1,1,2,2,2,3,3,3))
  #################### SNR simulations
  set.seed(123)
  sim <- as.vector(MIC_sim(alpha = 1, nsub = 1, Ct = c(3), fs = 200, segs = 1, iSNR=0)$Data[[1]])
  spec1 <- spec.pgram(sim,plot = F)$spec[1:50]; spec2 <- spec.parzen(sim,a = 100,nn = 100)$spec[1:50]
  plot(spec1/sum(spec1)); lines(spec2/sum(spec2), col = 'red')
  set.seed(123)
  sim <- as.vector(MIC_sim(alpha = 1, nsub = 1, Ct = c(3), fs = 200, segs = 1, iSNR=1)$Data[[1]])
  spec3 <- spec.pgram(sim,plot = F)$spec[1:50]; spec4 <- spec.parzen(sim,a = 100,nn = 100)$spec[1:50]
  points(spec3/sum(spec3), col = 'blue'); lines(spec4/sum(spec4), col = "blue")

  #################### To get the plotting coordinates
  coords <- read.csv("~/Dropbox/EEG_files/Data/EEG_coordinates.csv", header = T)
  coords <- coords[-129,]
  coords$theta <- coords$theta/180*pi;coords$phi <- coords$phi/180*pi

  coords <- coords[c(-125:-128),]

  off <- 16.5
  c <- sqrt(coords$radius^2+off^2-2*off*coords$radius*cos(pi/2+coords$phi))
  newr <- tan(acos((off^2+c^2-coords$radius^2)/(2*off*c)))
  newr[c(48,119)] <- .90*newr[c(48,119)]; newr <- newr/max(newr + .05)

  coords$x <- newr*sin(coords$theta); coords$y <- newr*cos(coords$theta)
  plot(coords$x, coords$y,xlim=c(-1,1),ylim = c(-1,1))
  # points(y=(newr*cos(coords$theta))[125:129],x=(newr*sin(coords$theta))[125:129], col = 'blue')
  EEGcoords <- coords[,c(1,5,6)]; devtools::use_data(EEGcoords, internal = T)
}
