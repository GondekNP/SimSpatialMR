sim.secr.replicate <- function(trapPath = 'detectorfileScaled.csv', wd = 'C:/Users/Nick/Documents/SubSECR/', baseLogsDir = 'ErrorLogs',
                           baseRdsDir = 'RepObjects', summaryPath = 'SimSECR_Summary.csv', parList, subTypes = c('full', 'SRS25'),
                           g0Model = c("g0 ~ 1")){
  #Function to automate the process of attempting a SECR run and saving the eventual error
  #or the successful output (NULL) in order to distribute to cores/instances
  setwd(wd)

  library(digest)
  library(sp)
  library(rgeos)
  library(adehabitatHR)
  library(tidyverse)
  library(doParallel)
  
  library(secr)
  library(secrdesign)
  gc()
  #initialize par conditions - in REAL PARAMETER terms 
  g0Par = parList$g0Par
  dPar = parList$dPar
  sigPar = parList$sigPar
  fitBuffPar = parList$fitBuffPar
  popnBuffPar = parList$popnBuffPar
  #establish where error logs will be saved if any 
  logsDir = file.path(baseLogsDir, paste0("dPar_", dPar),
                         paste0("g0Par_", g0Par),
                         paste0("sigPar_", sigPar),
                         paste0('buffPar_',fitBuffPar))
  if(!dir.exists(logsDir)){dir.create(logsDir, recursive = TRUE)}
  rdsDir = file.path(baseRdsDir, paste0("dPar_", dPar),
                      paste0("g0Par_", g0Par),
                      paste0("sigPar_", sigPar),
                      paste0('buffPar_',fitBuffPar),
                      paste0('model_',g0Model))
  if(!dir.exists(rdsDir)){dir.create(rdsDir, recursive = TRUE)}
  
  trapGrid <- read_csv(trapPath,
                       col_names = FALSE) %>% dplyr::select(X2, X3)
  colnames(trapGrid) <- c("x","y")
  trapGrid <- read.traps(data=trapGrid, detector = 'multi')
  
  #initialize lists to organize subsampled and full capthists and hairsamps data frames
  hairSamps <- list()
  capthists <- list()
  modFit <- list()

  trialPop <- sim.popn(D = dPar, core = trapGrid, buffer = popnBuffPar, model2d = "IHP")
  capthists$full <- sim.capthist(trapGrid, popn = trialPop, detectpar = list(g0 = qlogis(g0Par), 
                                                                             sigma = log(sigPar)), noccasions = 6)
  simID = digest(capthists$full) #cryptographic uniqueID based on capture history. Same capthist will have same hash 
  
  for (subType in subTypes){
    errPath = file.path(logsDir, paste0(simID, "_", subType, "_errorText.txt"))
    # browser()
    tryCatch({
          hairSamps$full <- as.data.frame(capthists$full) #full will be part of all reps
      
          #subsample (or dont), save as capthist and fit selected model
          if(!str_detect(subType, "full")){
            hairSamps[[subType]] <- hairSamps$full %>% sample_n(floor(nrow(hairSamps$full)*(as.numeric(str_extract(subType, "\\d\\d"))/100)))
          } #only robust to srs now
                                            
          capthists[[subType]] <- make.capthist(captures = hairSamps[[subType]], traps = trapGrid, fmt = 'trapID')
          modFit[[subType]] <- secr.fit(capthist = capthists[[subType]], buffer = fitBuffPar, detectfn = 0, model = as.formula(g0Model))
          
          #get derived parameters if successful
          g0Hat = plogis(modFit[[subType]]$fit$estimate[2])
          g0HatBias = g0Hat / g0Par
          
          sigHat = plogis(modFit[[subType]]$fit$estimate[3])
          sigHatBias = sigHat / sigPar
          
          dHat = derived(modFit[[subType]])["D","estimate"]
          dHatBiasPct = dHat/dPar
          dHatBias = dHat - dPar
          
          errorText = "SUCCESS"
          write_lines(errorText, path = errPath)
        }, error = function(e){
           write_lines(as.character(e), path = errPath)
        })
        errorText = paste(read_lines(errPath), collapse = '\n')#not sure how to avoid this step... looks like the error handling makes its own environment
        
        if (!is.null(hairSamps$full)) {
          nFull = nrow(hairSamps$full)
          nSub = nrow(hairSamps[[subType]])
          subFullProp = nSub/nFull
        }
        
        for (i in c('g0Hat', 'g0HatBias', 'sigHat', 'sigHatBias', 'dHat', 'dHatBias', 'dHatBiasPct', 'g0Hat', 'g0HatBias', 'nFull', 'nSub', 'subFullProp')){
          if (!exists(i)) {assign(i, NA)}
        }
        objPath = file.path(rdsDir, paste0(simID, '_', subType, '.rds'))
        #determine useful info from reps and save in data frame/list backup
        newLine = data.frame(simID, subType, dPar, dHat, dHatBiasPct, dHatBias, sigPar, 
                             sigHat, sigHatBias, g0Par, g0Hat, g0HatBias,
                             nFull, nSub, subFullProp, errorText, objPath)
        newObj = list(simID, subType, dPar, dHat, dHatBiasPct, dHatBias, sigPar, 
                            sigHat, sigHatBias, g0Par, g0Hat, g0HatBias,
                            nFull, nSub, subFullProp,
                            fullSamps = ifelse(is.null(hairSamps$full), NA, hairSamps$full),
                            subSamps = ifelse(is.null(hairSamps[[subType]]), NA, hairSamps[[subType]]),
                            fitted = ifelse(is.null(modFit[[subType]]), NA, modFit[[subType]]),
                            errorText, objPath)
        saveRDS(newObj, objPath)
        ifelse(file.exists(summaryPath), write_csv(newLine, summaryPath, append = T), write_csv(newLine, summaryPath, append = F))
  
    }
  
}