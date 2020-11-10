set.seed(52246)
library(magrittr)
'%!in%' <- function(x,y)!(x %in% y)
## Run Tests
setwd("C:/Users/Matthew/Documents/Courses/Kai/Final Results/Final Plots/DifferentialAbundanceTests")
dataWD <- "C:/Users/Matthew/Documents/Courses/Kai/Final Results/Final Plots/DataGeneration/"
source("Tests.R")
source("Normalizations.R")
source("RunTestsFunction.R")
source("fastTTest.R")


library(phyloseq)
library(ALDEx2)
library(limma)
library(edgeR)
library(DESeq2)
library(metagenomeSeq)
library(reshape2)
#library(doParallel)

  ## Load Data
AllData <- list()
load(paste0(dataWD,"BrokenPromiseData_NegBin.RData"))
AllData <- c(AllData,BrokenPromiseData)
load(paste0(dataWD,"BrokenPromiseData_BetaBin.RData"))
AllData <- c(AllData,BrokenPromiseData)

startt <- Sys.time()
MasterRes <- data.frame()
pvalList <- list()

seq <- 1:length(AllData)
MasterRes <- data.frame()

## Run Tests 
for(dat in seq){
  
  print(dat)
  datalist <- AllData[[dat]]
  physeq <- datalist$physeq
  truede <- datalist$truede
  degenes <- rownames(otu_table(physeq))[truede]
  
  if(length(degenes) == 0) degenes <- c()
  sample_data <- as(sample_data(physeq), "data.frame")
  if(is.null(sample_data$condition)){
    sample_data$condition <- sample_data$group
  } else if(is.null(sample_data$group)){
    sample_data$group <- as.character(sample_data$condition)
  }
  sample_data$group <- as.character(sample_data$group)
  sample_data$condition <- as.character(sample_data$condition)
  sample_data(physeq) <- sample_data
  
  ## returns Test, Norm Sens, Spec, FDR
  ## this is a matrix 
  Res <- try(runTests(physeq, truede = truede, degenes = degenes))
  if(class(Res) != "try-error"){
    ResTemp <- data.frame("test" = NULL, "norm" = NULL, "sens" = NULL, "spec" = NULL, "fdr" = NULL, stringsAsFactors = FALSE)
    for(resentry in 1:(length(Res))){
      ResEntry <- Res[[resentry]]$res
      ResTemp <- rbind(ResTemp, ResEntry)
    }
    
    ## for the pvals 
    for(resentry in 1:(length(Res))){
      testname <- names(Res)[resentry]
      itername <- names(AllData)[dat]
      pv_name <- paste(testname, itername, sep = "/")
      pvalList[[pv_name]] <- unlist(Res[[resentry]]$pvals)
    }
    ResTemp$fdr <- sapply(ResTemp$fdr, function(x){
      if(x == "NaN"){
        0
      } else{
        as.numeric(as.character(x))
      }
    })
    ResTemp$sens <- sapply(ResTemp$sens, function(x){
      if(x == "NaN"){
        0
      } else{
        as.numeric(as.character(x))
      }
    })
      ## fix specificity 
    ResTemp$spec <- sapply(ResTemp$spec, function(x){
      if(x == "NaN"){
        0
      } else{
        as.numeric(as.character(x))
      }
    })
    Res <- ResTemp
    Res$norm <- tolower(as.character(Res$norm)) ## fixes norm bug
  
    ##
    nrow <- nrow(Res)
    newname <- (names(AllData[dat]))
    parameters <- strsplit(newname, "_")
   #sim <- rep(parameters[[1]][1], nrow)
    sim <- parameters[[1]][1]
    effect <- (parameters[[1]][2]) 
    m <- (parameters[[1]][3]) %>% trimws()
    m <- ncol(as(otu_table(physeq),"matrix"))/2
    propde <- (parameters[[1]][4])
    iter <- (parameters[[1]][5])
    sampleType <- parameters[[1]][6]
    ## 11/9/2019
      distr <- parameters[[1]][7]
      sim <- c(ifelse(is.na(distr), sim, 
                    ifelse(distr == "dirmult", "dirmult", sim)))
    ##  
    params <- (paste(Res[,"test"], Res[,"norm"], sim, effect, m, propde, sampleType, sep = "_"))
    
    ##
    Res <- cbind(Res, 
                 sim, effect, m, propde, iter, params)
    #Res
    MasterRes <- rbind(MasterRes, Res)
  } else{
    cat('\n\n\tWARNING: MISSED DATASET ON ', names(AllData)[dat],'\n\n')
  }
}
endd <- Sys.time()
endd - startt
MasterRes$spec <- as.numeric(as.character(MasterRes$spec))
MasterRes$fdr <- as.numeric(as.character(MasterRes$fdr))
MasterRes$m <- as.numeric(as.character(MasterRes$m))
MasterRes$effect <- as.numeric(as.character(MasterRes$effect))
MasterRes$propde <- as.numeric(as.character(MasterRes$propde))
MasterRes$test <- as.character(MasterRes$test)
MasterRes$norm <- as.character(MasterRes$norm)
MasterRes$sim <- as.character(MasterRes$sim)
MasterRes$test_norm <- paste(MasterRes$test, MasterRes$norm, sep = "_")
save(MasterRes,file = paste0("MasterRes",".RData"))
save(pvalList,file = paste0("pvalList",".RData"))




