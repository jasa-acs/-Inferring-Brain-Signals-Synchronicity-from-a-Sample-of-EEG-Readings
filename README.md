# Inferring-Brain-Signals-Synchronicity-from-a-Sample-of-EEG-Readings

# Author Contributions Checklist Form

## Data

### Abstract 

The data consist of Electroencephalography (EEG) recordings from participants who were asked to watch videos of non-social images on a computer monitor for 2 to 6 minutes. The recording was simultaneous through multiple (129) electrodes, at a sampling frequency of 250Hz, and subsequently segmented into units of 1024ms which correspond to 256 observations per segment. Therefore, the data structure for each participant is a 3D-array of dimensions 129 × 256 × no. of segments, where the number of available segments vary by participants.

### Availability 

The data in the case study is owned by Dr. Jeste’s group at UCLA. Our collaborators are currently working on several other manuscripts involving the same study. As the study itself represents the major scientific capital for the group, we do not plan on making the data publicly available at the moment.

## Code

### Abstract

The MIC2 (Multilevel Integrative Clustering) package performs integrative clustering on highly structured data. Our development includes functions for the simulation of time series, MCMC simulation from the target posterior, as well as functions implementing model search and post- processing of posterior samples.

The MIC2 package was developed in the R statistical programming framework and distributed under GPL-2.
Version 2.1.1 is freely available for download at [https://github.com/Qian-Li/MIC2](https://github.com/Qian-Li/MIC2.).

### Description

Dependencies include Rcpp, RcppArdmadillo, Class.

## Instructions for Use

### Reproducibility 

The data analysis workflow is summarized in the file Readme.md as part of the MIC2 documentation on GitHub. More details are available in the package vignette. All simulation results are easily reproduced using the function MIC_sim and related documentation in the package.

### Replication 

All functions in the MIC2 package are documented with usage examples.
