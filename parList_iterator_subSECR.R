library(secr)
library(secrdesign)
library(digest)
library(sp)
library(rgeos)
library(adehabitatHR)
library(tidyverse)
library(doParallel)


source('~/SubSECR/robust_replicate_function_subSECR.R')


cl = makeCluster(detectCores()-3)
registerDoParallel(cl)

foreach(i = 1:600) %dopar% {
  dPars = seq(525.25, 5252.5, length.out = 10)
  g0Pars = seq(.5, .8, length.out = 10)
  sigPars = seq(10,1,-1)
  parLists <- list()
  j = 1
  for (k in 1:10){
    for (l in 1:10){
      for (o in 1:10){
        parLists[[j]] <- list(g0Par = g0Pars[k], dPar = dPars[l], sigPar = sigPars[o], popnBuffPar = 30, fitBuffPar = 10) 
        j = j+1
      }
    }
  }
  
  for (p in 1:10000){
    for (k in 1:length(parLists)){
      sim.secr.replicate(parList = parLists[[k]])
    }
  }
}