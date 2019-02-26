# MIC2
Multilevel Integrative Clustering (MIC recoded in C++), packaged by Qian Li

# Installation:
```r
devtools::install_github("Qian-Li/MIC2")
```

More technical details are available in [this introduction](MIC2-intro.md) and our [manuscript](https://arxiv.org/abs/1609.09532).

# Quick Example

Time series simulation:
```r
library(MIC2)
## 3 subjects, each has 10 segments of simulated signals, coherence = 0.9
ts_sim <- MIC_sim(alpha = 0.9, nsub = 3, segs = 10, fs = 200)
```

Data preparation:
```r
## Spectral Estimation and Eigen-Laplacian transformation
list_data <- lapply(ts_sim$Data, function(x) MIC_prep(X = x, d = 4,
                                  par.spec = c(80, 50, 100), par.win = c(3, 1)))
```


MIC clustering:
```r
## MIC fitting with 10k iterations (9 secs)
start <- Sys.time()
output <- MIC(data = list_data, K = 4, nit = 10000); 
Sys.time() - start
```

Clustering result:
```r
## Group label accuracy:
sum(align(ts_sim$C, output$clustS, type = 'vec') == ts_sim$C) / 40

## Individual labels accuracy:
true_c <- lapply(1:3,function(i) ts_sim$Ci[,i]); est_c <- lapply(1:3,function(i) output$clustC[i,])
mapply(function(x,y) sum(align(x,y,'vec')==x)/40,
    true_c, est_c)
```
