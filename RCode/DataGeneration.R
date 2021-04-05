

## I'm using R version 3.5.1 right now

## Contains code necessary to generate simulated datasets considered
## Code was adapted for non-parallel computation and streamlined for H1 midVagina Template
## Source code from https://users.ugent.be/~shawinke/ABrokenPromise/02_dataGeneration.html  

##########################################################################################
set.seed(52246)
dataWD <-
"C:/Users/Matthew/Documents/Courses/Kai/Final Results/Final Plots/DataGeneration"
setwd(dataWD)
source("msWaldHMP.R") ## copy/pasted from github
knitr::opts_chunk$set(
  cache = FALSE,
  tidy = TRUE,
  autodep = FALSE,
  root.dir = WD,
  eval = TRUE
)

# The required package list:
reqpkg = c("parallel",
           "phyloseq",
           "MASS",
           "SpiecEasi",
           "magrittr",
           "TailRank")
# Load all required packages and show version
for (i in reqpkg)
{
  print(i)
  print(packageVersion(i))
  library(
    i,
    quietly = TRUE,
    verbose = FALSE,
    warn.conflicts = FALSE,
    character.only = TRUE
  )
}


# # distribution to generate data from
#distribs = c("betabinCor")
distribs = c("negbinCorOut")

# # number of repeat datasets per unique combo of parameters
reps <- 1:25L

# # true positive rate
TPR <- .1

# # for labelling the file later
TPR_label <- "1"

# # for labelling the file later
letter <- "A"

# # Minimum number of reads to consider an OTU 'observed' in a sample
minReads <- 1L

# # The delimiter in the command parameter string
delim <- "_"

# # Define the different biological source templates to use
sampleTypes <- c("Mid.vagina")

# # Define the ceiling in the number of OTUs to consider in the template

nOTUs <- 1000L

# # Define the number of samples in each class of a simulated experiment
nObs <- c(5,15,25)
# nObs <- c(5, 100)

# # The different values of effec5t size to apply
foldEffect <- c(3,5)

# # The number of cores used in parallel computing (I didn't use any)
nCores <- 1

# # The covariance estimation method
covEstMethod = "glasso"

# # Biologically relevant variables
variables = c("IBDbin", "Penbin", "Sexbin")
plasmSampleNames = c("Stool_sex", "Tongue.dorsum_sex")
nObs = sort(nObs, decreasing = TRUE)
# # Define the simulation parameters combinations
simParams <-
  apply(expand.grid(sampleTypes, reps, nObs, distribs),
        1L,
        paste,
        collapse = delim)
simParams <-
  gsub(
    pattern = " ",
    replacement = "",
    x = simParams,
    fixed = TRUE
  )
simParamsLabels <-
  c("SampleType", "Replicate", "nSamples", "Distribution")
simParamsH1 <-
  apply(expand.grid(sampleTypes, reps, foldEffect, nObs, distribs),
        1L,
        paste,
        collapse = delim)
simParamsH1 <-
  gsub(
    pattern = " ",
    replacement = "",
    x = simParamsH1,
    fixed = TRUE
  )

# Define the labels to go with each element of the simulation parameter
# after splitting on the delimiter
simParamsLabelsH1 <-
  c("SampleType",
    "Replicate",
    "EffectSize",
    "nSamples",
    "Distribution")

seq <- 1:length(simParamsH1)
#####################################################

load("physeqListV13.RData")
load("physeqListV35.RData")
physeqListV13AG <- physeqListV13
if (!file.exists("physeqList4Trim.RData")) {
  OTUsKeep = lapply(physeqListV13AG, function(x) {
    relAbundances = taxa_sums(x)
    names(sort(relAbundances, decreasing = TRUE)[1:nOTUs])
  })
  physeqList4Trim = mapply(
    physeqListV13AG,
    OTUsKeep,
    FUN = function(phy, otu) {
      prune_taxa(phy, taxa = otu)
    }
  )
  #physeqList4Trim$AGstool <- NULL
  #rm(physeqListV13AG)
  save(physeqList4Trim, file = "physeqList4Trim.RData")
} else {
  load("physeqList4Trim.RData")
}
print("Done Loading Data")

if (!file.exists(file = "piMoMs.RData")) {
  piMoMs <- lapply(physeqList4Trim, function(x) {
    if (taxa_are_rows(x)) {
      piMoM4Wald(t(x@otu_table@.Data))
    } else {
      piMoM4Wald(x@otu_table@.Data)
    }
  })
  thetaMoMs <- sapply(physeqList4Trim, function(x) {
    if (taxa_are_rows(x)) {
      weirMoM4Wald(t(x@otu_table@.Data), se = FALSE)
    } else {
      weirMoM4Wald(x@otu_table@.Data)
    }
  })
  save(thetaMoMs, piMoMs, file = "piMoMs.RData")
} else {
  load(file = "piMoMs.RData")
}

#### Negative binomial parameter estimation ####
if (!file.exists("MLES.RData")) {
  #clu <- makeCluster(nCores, outfile = "logFileNBfits.txt")
  #clusterEvalQ(cl = clu, {
  require(phyloseq, quietly = TRUE)
  require(MASS, quietly = TRUE)
  #})
  
  NBfitsList <- lapply(physeqList4Trim, function(x) {
    if (taxa_are_rows(x)) {
      logLibSizes = log(colSums(x@otu_table@.Data))
      apply(x@otu_table@.Data, 1, function(y) {
        try(glm.nb(y ~ offset(logLibSizes), link = "log"), silent = TRUE)
      })
    } else {
      logLibSizes = log(rowSums(x@otu_table@.Data))
      apply(x@otu_table@.Data, 2, function(y) {
        try(glm.nb(y ~ offset(logLibSizes), link = "log"), silent = TRUE)
      })
    }
  })
  #stopCluster(clu)
  rhoMLEs = lapply(NBfitsList, function(x) {
    tmp = sapply(x, function(y) {
      if (class(y)[1] != "negbin") {
        NA
      } else {
        exp(y$coef[1])
      }
    })
    names(tmp) = names(x)
    res = tmp[!is.na(tmp)]
    res / sum(res)
  })  #Renormalize!
  
  
  phiMLEs = lapply(NBfitsList, function(x) {
    tmp = sapply(x, function(y) {
      if (class(y)[1] != "negbin") {
        NA
      } else {
        1 / y$theta
      }
    })
    names(tmp) = names(x)
    tmp[!is.na(tmp)]
  })
  PearRes = lapply(NBfitsList, function(x) {
    pr <- try(sapply(x, residuals, type = "pearson"), silent = TRUE)
    if (class(pr) == "try-error") {
      invisible()
    } else{
      return(pr)
    }
  })
  save(list = c("rhoMLEs", "phiMLEs", "PearRes"),
       file = "MLES.RData")
} else {
  load("MLES.RData")
}

ExtrNBouts = function(PearRes, PearsonCutOff = 5) {
  outliers = abs(PearRes) > PearsonCutOff
  freqVec = rowSums(outliers) / ncol(outliers)  #Relative frequency: outliers per taxon
  PearVec = PearRes[outliers]
  list(freqOut = freqVec, Pres = PearVec)
}

## Get pearson residuals for introducing outliers later
PearRes <- PearRes#[-(which(sapply(PearRes,is.null),arr.ind=TRUE))]
OutLieList = lapply(PearRes, ExtrNBouts)
save(OutLieList, file = "outLieList.RData")
print("Done with Outliers")
print(length(OutLieList))

# if (!file.exists("CovListEst.RData")) {
#   covListEst = lapply(physeqList4Trim, spiec.easi, icov.select.params = list(ncores = nCores))
#   save(covListEst, file = "CovListEst.RData")
# } else {

## Note From Matt: Running the actual code broke two of my university's computers
## I requested this directly from Hawinkel et al, who kindly provided it instead.
load(file = "covList.RData")

### generate Dirichlet realisations, taken from gtools (identical in MCMCpack)
rDirichlet <- function(n, alpha) {
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  res <- x / as.vector(sm)
  res[res <= 10 * .Machine$double.eps] <- 0
  res
}

# # A custom quantile beta-binomial function with `na.rm=TRUE`. Still relies
# on the Tailrank package
qbetabin = function(p, N, u, v) {
  pp <- cumsum(dbb(0:N, N, u, v))
  sapply(p, function(x)
    sum(pp < x, na.rm = TRUE))
}

# # A function to generate correlated multivariate betabinomial data pi: a
# vector of proportions, summing to 1 libSizes: library sizes theta: the
# overdispersion parameter
rmvbetabin = function(n, pi, Sigma, theta, libSizes, ...) {
  Sigma <- as.matrix(Sigma)
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  if (missing(pi))
    stop("pi is required")
  if (length(pi) != dim(Sigma)[1])
    stop("Sigma and pi dimensions don't match")
  if (missing(theta)) {
    stop("No overdispersion parameter supplied")
  }
  d <- length(pi)
  normd <-
    rmvnorm(n, rep(0, d), Sigma = Cor)  #The normal-to-anything framework
  unif <- pnorm(normd)
  data <-
    mapply(
      unif,
      rep(pi, each = nrow(unif)),
      libSizes,
      FUN = function(u,
                     p, l) {
        alphaPar = p * (1 - theta) / theta
        betaPar = (1 - p) * (1 - theta) / theta
        qbetabin(u, N = l, u = alphaPar, v = betaPar)
      }
    )
  data <- .fixInf(data)
  return(data)
}

# First an auxiliary function
.fixInf <- function(data) {
  # hacky way of replacing infinite values with the col max + 1
  if (any(is.infinite(data))) {
    data <- apply(data, 2, function(x) {
      if (any(is.infinite(x))) {
        x[ind <- which(is.infinite(x))] <- NA
        x[ind] <- max(x, na.rm = TRUE) + 1
      }
      x
    })
  }
  data
}
# # Generate correlated NB data, given a covariance matrix n: number of
# observations mu: means of NB distribution Sigma: a positive definite
# covariance matrix ks: overdispersion parameters (size)
rmvnegbin = function(n, mu, Sigma, ks, ...) {
  Sigma <- as.matrix(Sigma)
  Cor <- cov2cor(Sigma)
  SDs <- sqrt(diag(Sigma))
  if (missing(mu))
    stop("mu is required")
  if (dim(mu)[2] != dim(Sigma)[2])
    stop("Sigma and mu dimensions don't match")
  if (missing(ks)) {
    ks <-
      unlist(lapply(1:length(SDs), function(i)
        .negbin_getK(mu[i], SDs[i])))
  }
  d <- dim(mu)[2]
  normd <-
    rmvnorm(n, rep(0, d), Sigma = Cor)  #The normal-to-anything framework
  unif <- pnorm(normd)
  data <- t(qnbinom(t(unif), mu = t(mu), size = ks))
  data <- .fixInf(data)
  return(data)
}


### `sampleSizes` is a vectors, `alphas` ,`phis`, 'rhos' and `Sigma` matrices, `libSizes` a list
### final matrix has as rownames the sample names taken from `libSizes`
### and as colnames OTU names taken from rownames of `alphas` or `rhos`
### distribution is the
countsGen <-
  function(sampleSizes,
           distribution = c(
             "negbinNoCor",
             "negbinCor",
             "dirmult",
             "betabinCor",
             "negbinCorOut",
             "negbinNoCorOut",
             "betabinCorOut"
           ),
           alphas = NULL,
           theta = NULL,
           rhos = NULL,
           phis = NULL,
           libSizes = NULL,
           Sigma = NULL,
           onlyCounts = TRUE,
           outLiers = NULL)
  {
    if (!is.list(libSizes))
    {
      stop("`libSizes` must be a list of length `length(sampleSizes)`")
    } else {
    }
    
    libSizes <- unlist(libSizes, use.names = TRUE)
    if (distribution %in% c("negbinCorOut", "negbinNoCorOut", "betabinCorOut") &
        is.null(outLiers))
    {
      stop("No outlier matrix supplied")
    }
    if (!distribution %in% c(
      "negbinNoCor",
      "negbinCor",
      "dirmult",
      "betabinCor",
      "negbinCorOut",
      "negbinNoCorOut",
      "betabinCorOut",
      "betabinCorOut"
    ))
    {
      stop("No valid count distribution supplied")
    } else if (distribution %in% c("negbinCor",
                                   "negbinNoCor",
                                   "negbinCorOut",
                                   "negbinNoCorOut",
                                   "betabinCorOut"))
      ## Negative binomial
    {
      if (is.null(rhos) | is.null(phis))
      {
        stop("No valid NB parameters supplied")
      } else{
      }
      
      nbData <-
        matrix(NA_integer_,
               nrow = sum(sampleSizes),
               ncol = nrow(rhos)) #All datasets have the same number of OTUs
      samNames <-
        rep(paste0("grp", seq_along(sampleSizes)), sampleSizes)
      samNames <-
        paste(samNames, rep.int(seq_len(sampleSizes[1]), length(sampleSizes)),
              sep = ":")
      
      rownames(nbData) <- samNames
      colnames(nbData) <- rownames(as.matrix(rhos))
      samSizeSeq <- c(0L, cumsum(sampleSizes))
      if (distribution %in% c("negbinNoCor", "negbinNoCorOut"))
      {
        for (nRun in seq_along(sampleSizes))
        {
          ## selected indices to generate
          indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
          
          ## Negative binomial draws
          nbData[indSel,] <-
            mapply(
              rhos[, nRun],
              phis[, nRun],
              FUN = function(rho, phi)
              {
                rnbinom(n = sampleSizes[nRun],
                        mu = rho * libSizes[indSel],
                        size = 1 / phi)
              }
            )
        }
      } else if (distribution %in% c("negbinCor", "negbinCorOut"))
      {
        if (is.null(Sigma))
        {
          stop("No correlation matrix given")
        } else if (dim(Sigma)[1] != dim(Sigma)[2])
        {
          stop("Correlation matrix is not square")
        } else{
        }
        
        for (nRun in seq_along(sampleSizes))
        {
          ## selected indices to generate
          indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
          
          ## Negative binomial draws with underlying correlation
          
          nbData[indSel,] <-
            rmvnegbin(
              n = sampleSizes[nRun],
              mu = t(tcrossprod(rhos[, nRun], libSizes[indSel])),
              ks = 1 / phis[, nRun],
              Sigma = Sigma
            )
        } #end: for
      } #end- if negbinCor
      if (distribution %in% c("negbinCorOut", "negbinNoCorOut", "betabinCorOut")) {
        #Introduce outliers into generated data
        
        #Introduce outliers randomly over entire counts matrix
        nSamples = dim(nbData)[1]
        nTaxa = dim(nbData)[2]
        outFracs = outLiers[["freqOut"]][names(libSizes)[names(libSizes) %in% names(outLiers[["freqOut"]])]] #Fraction of outliers in each sample, kept connected with libSizes
        nOuts = rbinom(nSamples, nTaxa, outFracs) #Number of outliers in each sample
        nbData = t(sapply(1:nSamples,  function(i) {
          if (nOuts[i] == 0) {
            return(nbData[i, ])
          }
          pearRes = sample(outLiers[["Pres"]], nOuts[i], replace = TRUE) #Sample Pearson residuals
          taxaIDs = sample(colnames(nbData), nOuts[i], replace = FALSE)
          expects = libSizes[i] * rhos[taxaIDs, 1] #Expected outcomes
          newValues = sapply(round(sqrt(expects * (
            1 + expects * phis[i, 1]
          )) * (pearRes) + expects), max, 0) #Reconstruct outliers from Pearson residuals. Round and set negative values to zero
          nbData[i, taxaIDs] = newValues
          nbData[i, ]
        }))
        rownames(nbData) <- samNames
      } else{
      }
      nbData
    } else if (distribution %in% c("dirmult", "betabinCor")) {
      if (length(sampleSizes) != NCOL(alphas))
      {
        stop("length(sampleSizes) must be the same of ncol(alphas)")
      } else {
      }
      
      dmData <-
        matrix(NA_integer_,
               nrow = sum(sampleSizes),
               ncol  = dim(Sigma)[1])
      samNames <-
        rep(paste0("grp", seq_along(sampleSizes)), sampleSizes)
      #   samNames <- paste(samNames, names(libSizes), sep = ":")
      samNames <-
        paste(samNames, rep.int(seq_len(sampleSizes[1]), length(sampleSizes)),
              sep = ":")
      #   Ntaxa=dim(Sigma)[1]
      #   id=sample(1:nrow(alphas), Ntaxa)
      #   alphas=alphas[id,]
      alphas = alphas / sum(alphas[, 1]) #renormalize based on the unchanged alphas
      rownames(dmData) <- samNames
      colnames(dmData) <- rownames(alphas)
      piDir <- dmData
      piDir[] <- NA_real_
      samSizeSeq <- c(0L, cumsum(sampleSizes))
      
      if (distribution == "dirmult")
      {
        ## gamma parameter for Dirichlet distribution
        gammas <- alphas * (1 - theta) / theta
        for (nRun in seq_along(sampleSizes))
        {
          ## selected indices to generate
          indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
          
          ## Dirichlet draws
          piDir[indSel,] <- rDirichlet(# piDir[indSel, ] <- gtools:::rdirichlet(
            n = sampleSizes[nRun], alpha = gammas[, nRun])
          
          ## Multinomial draws with Dirichlet probabilities (Dirichlet-Multinomial)
          dmData[indSel,] <- t(sapply(indSel,
                                      function(iRun)
                                      {
                                        rmultinom(n = 1L,
                                                  size = libSizes[iRun],
                                                  prob = piDir[iRun,])
                                      }))
        }# END - for: loop along sample sizes
      } else if (distribution == "betabinCor")
      {
        if (is.null(Sigma))
        {
          stop("No correlation matrix given")
        } else if (dim(Sigma)[1] != dim(Sigma)[2]) {
          stop("Correlation matrix is not square")
        } else{
        }
        for (nRun in seq_along(sampleSizes))
        {
          ## selected indices to generate
          indSel <- (samSizeSeq[nRun] + 1L):samSizeSeq[nRun + 1L]
          
          dmData[indSel,] <-
            rmvbetabin(
              n = sampleSizes[nRun],
              pi = alphas[, nRun],
              libSizes = libSizes[indSel],
              Sigma = Sigma,
              theta = theta
            )
        }
      }
      
      if (onlyCounts)
      {
        dmData
      } else
      {
        list("dmData" = dmData, "piDir" = piDir)
      }
    }# END - if: distributions
  }# END - function: countsGen


#minReads = 1 #- an OTU is considered present in a sample if it has at least 1 read
#prevalence = 0.05 #- an OTU is kept if it is present in at least 5% of samples
#prevalence = 0
# # Trim by prevalence and total OTU reads
# ** ** NOTE: MATT SET MIN READS TO 0 ** **
simpleTrimGen <- function(obj,
                          minReads = 1L,
                          minPrev = .05) {
  # `prevalence` is the fraction of samples in which an OTU is observed at
  # least `minReads` times.
  if (class(obj) == "phyloseq") {
    taxRows <- taxa_are_rows(obj)
    if (!taxRows) {
      obj <- t(obj)
    } else {
      
    }
    otuTab <- as(otu_table(obj), "matrix")
  } else {
    otuTab <- obj
  }  # END - ifelse: obj is *phyloseq* or just *matrix*
  
  ## sort OTUs first by prevalence, and then by total reads per OTU
  prevalence <- rowMeans(otuTab >= minReads)
  #prevalence <- rowMeans(otuTab > minReads)
  
  ## Matt changed this: 
  ## Since ALDEx2 and rank normalization need full datasets, 
  ## we will filter 5% later for only the other methodologies, and keep full 
  ## datasets here
  indOTUs2Keep <- (prevalence > 0)
  #indOTUs2Keep <- (prevalence >= minPrev)
  
  if (class(obj) == "phyloseq") {
    obj = prune_taxa(obj, taxa = indOTUs2Keep)
    return(obj)
  } else {
    return(otuTab[indOTUs2Keep,])
  }
}  # END - function: simpleTrim general

# # Check if less than _minOTUs_ taxa are present in each sample

fewTaxa <- function(physeq, minOTUs = 0L) {
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  } else {
    
  }
  
  any(colSums(otu_table(physeq) > 0, na.rm = TRUE) < minOTUs)
}


addFoldChange = function(rhos,
                         fc,
                         H1frac = TPR,
                         compensate = FALSE) {
  if (fc == 1) {
    return(rhos)
  }
  nTaxa = length(rhos)
  if (compensate) {
    nOTUsUp = round(nTaxa * H1frac * (1 / (fc + 1)))  #Upregulated taxa
    nOTUsDown = round(nTaxa * H1frac - nOTUsUp)  #Downregulated taxa
    # cond=TRUE while(cond){
    OTUids = sample(names(rhos), nOTUsUp + nOTUsDown, replace = FALSE)
    OTUidUps = OTUids[1:nOTUsUp]
    OTUidDowns = OTUids[(nOTUsUp + 1):(nOTUsDown + nOTUsUp)]
    rhos[OTUidUps] = rhos[OTUidUps] * fc  # Add fold change up
    rhos[OTUidDowns] = rhos[OTUidDowns] * (1 - sum(rhos[OTUidUps]) - sum(rhos[!(names(rhos) %in%
                                                                                  OTUids)])) /
      sum(rhos[OTUidDowns])  #And compensate the downs. This way the average FC is 5 in both directions and the TN taxa are really left untouched
    indTPup <- names(rhos) %in% OTUidUps
    newTaxaNamesUp <- paste0(names(rhos)[indTPup], "-TPup")
    indTPdown <- names(rhos) %in% OTUidDowns
    newTaxaNamesDown <- paste0(names(rhos)[indTPdown], "-TPdown")
    names(rhos)[indTPup] <- newTaxaNamesUp
    names(rhos)[indTPdown] <- newTaxaNamesDown
  } else {
    nOTUs = round(nTaxa * H1frac)  #DA taxa
    OTUids = sample(names(rhos), nOTUs, replace = FALSE)
    rhos[OTUids] = rhos[OTUids] * fc  # Add fold change up
    indTP <- names(rhos) %in% OTUids
    newTaxaNames <- paste0(names(rhos)[indTP], "-TPup")
    names(rhos)[indTP] <- newTaxaNames
  }
  rhos / sum(rhos)  #Renormalize.
}


microbioSim <-
  function(postfix,
           template,
           estPi,
           estTheta,
           nObs,
           estPhis,
           Covar,
           distrib,
           estRhos,
           outLiers,
           foldChange = 1,
           compensate = FALSE) {
    # Generate `nObs` simulated microbiomes with `libSizes` total reads each
    # where `libSizes` is a list where each element contains a vector of length
    # equal to the value of the corresponding element of `nObs`.  `libSizes`
    # contains samples drawn from `template` total reads.  `postfix` is a dummy
    # idenitifer added to help distinguish simulated samples in downstream code.
    libSizesOrig <- as(sample_sums(template), "integer")
    libSizes <-
      list(
        sample(libSizesOrig, size = nObs, replace = TRUE),
        sample(libSizesOrig,
               size = nObs, replace = TRUE)
      )
    
    # Actually create the simulated abundance table, both groups at once
    AltRhos = addFoldChange(estRhos, foldChange, compensate = compensate)  #Relative abundances of the other group
    defRhos = cbind(estRhos, AltRhos)
    rownames(defRhos) = names(AltRhos)
    AltAlphas = addFoldChange(estPi, foldChange)  #Relative abundances of the other group
    defAlphas = cbind(estPi, AltAlphas)
    rownames(defAlphas) = names(AltAlphas)
    
    counts <-
      countsGen(
        sampleSizes = c(nObs, nObs),
        alphas = defAlphas,
        theta = estTheta,
        onlyCounts = TRUE,
        libSizes = libSizes,
        rhos = defRhos,
        distribution = distrib,
        Sigma = Covar,
        phis = cbind(estPhis, estPhis),
        outLiers = outLiers
      )
    ## Add the OTU names to the OTU (column) indices, not needed with countsGen
    ## colnames(counts) <- taxa_names(template) Add new simulated sample_names to
    ## the row (sample) indices, not needed rownames(counts) <-
    ## paste(rownames(counts), '::', postfix, sep = '')
    
    # Put simulated abundances together with metadata as a phyloseq object, taxa
    # are rows here as it is consistent with other packages
    otuTab <- otu_table(t(counts), taxa_are_rows = TRUE)
    # Define data.frame that will become sample_data
    samNames <- sample_names(otuTab)
    samNames <-
      matrix(unlist(strsplit(
        x = samNames, split = ":", fixed = TRUE
      )),
      nrow = 2L)
    
    samData <-
      data.frame(
        group = samNames[1L,],
        sample = samNames[2L,],
        postfix = postfix,
        stringsAsFactors = FALSE
      )
    rownames(samData) <- sample_names(otuTab)
    samData <- sample_data(samData)
    # Return a phyloseq object
    return(phyloseq(otuTab, samData))
  }  # END - function: microbioSim


makePlasmodes = function(Sample,
                         samSize,
                         postfix,
                         var,
                         physeq,
                         meanLFDR = meanLFDR) {
  if (!taxa_are_rows(physeq)) {
    physeq = t(physeq)
  }
  counts = physeq@otu_table@.Data
  treatment = ifelse(sample_data(physeq)[[var]] == unique(sample_data(physeq)[[var]])[1],
                     "grp1",
                     "grp2")
  idNA = is.na(treatment)
  treatment = treatment[!idNA]
  counts = counts[,!idNA]
  if (max(table(treatment)) < 2 * samSize) {
    return(NULL)
  }
  # switchTrt = table(treatment)[1] < table(treatment)[2]
  plasmode = SimData(
    counts,
    treatment,
    sort.method = "unpaired",
    k.ind = samSize,
    n.diff = round(nrow(counts) * TPR),
    weights = meanLFDR[[paste(Sample,
                              var, sep = delim)]],
    n.genes = nrow(counts),
    norm.factors = colSums(counts)
  )
  # Construct phyloseq object from here
  otuTab = otu_table(plasmode$counts, taxa_are_rows = TRUE)
  samTab = sample_data(data.frame(
    group = ifelse(
      plasmode$treatment == unique(plasmode$treatment)[1],
      "grp1",
      "grp2"
    ),
    postfix = postfix
  ))
  physeq = phyloseq(samTab, otuTab)
  taxaNames = taxa_names(physeq)
  taxaNames[plasmode$DE.ind] = sapply(taxaNames[plasmode$DE.ind], paste0,
                                      "-TP")
  taxa_names(physeq) = taxaNames
  physeq
}
physeqListV13AG = lapply(sampleTypes, function(x) {
  physeqListV13AG[[x]]
})
names(physeqListV13AG) = sampleTypes

## Do the same for the HMP data based on gender and runcenter WUGC:
## Washington University Genome Center
HMPsubset = lapply(physeqListV13AG[names(physeqListV13AG) != "AGstool"], function(x) {
  prune_samples(x,
                samples = sample_data(x)$RUNCENTER == "WUGC" &
                  sample_data(x)$sex ==
                  "female")
})

cleanList0 = c(HMPsubset)

cleanList = lapply(cleanList0, function(x) {
  if (taxa_are_rows(x)) {
    piMoMs = piMoM4Wald(t(x@otu_table@.Data))
  } else {
    piMoMs = piMoM4Wald(x@otu_table@.Data)
  }
  names(piMoMs) = taxa_names(x)
  taxaKeep = names(sort(piMoMs, decreasing = TRUE)[1:nOTUs])
  prune_taxa(x = x, taxaKeep)
})
rm(cleanList0, HMPsubset)

subSample = function(physeq,
                     split = 0.5,
                     nObs,
                     replace = FALSE,
                     postfix) {
  if (nsamples(physeq) < 2 * nObs)
  {
    return(NULL)
  }  #If not enough samples to perform subsampling, return NULL
  idSample = sample(1:nsamples(physeq), nObs * 2, replace = replace)
  
  # A complication: we want to stick to the phyloseq framework, but then
  # samples cannot have identical names. We need a workaround when working
  # with sampling with replacement
  if (replace) {
    Table = table(idSample)
    if (max(Table) > 1) {
      for (x in (2:max(Table))) {
        samNumbers = as.integer(names(Table[Table >= x]))
        tmp = prune_samples(physeq, samples = sample_names(physeq)[samNumbers])
        # Now change the sample names to avoid conflicts
        sample_names(tmp) = paste0(sample_names(tmp), x)
        physeqTmp = merge_phyloseq(tmp, physeqTmp)
      }
    }
  } else {
    physeqTmp = prune_samples(physeq, samples = sample_names(physeq)[idSample])
  }
  groupSizes = round(c(split * nObs * 2, (1 - split) * nObs * 2))
  # Assign groups at random
  sample_data(physeqTmp) = data.frame(
    group = paste0("grp", as.integer(sample(c(
      rep(1,
          groupSizes[1]), rep(2, groupSizes[2])
    )))),
    sample = rep(seq(1, nObs),
                 2),
    postfix = postfix,
    row.names = sample_names(physeqTmp)
  )
  return(physeqTmp)
}


splitSample = function(physeq,
                       nEval,
                       postfix,
                       variable,
                       maxVerif = 100) {
  physeq = prune_samples(x = physeq, as.vector(!is.na(sample_data(physeq)[,
                                                                          variable])))  #Remove NA's
  samDataDF = data.frame(sample_data(physeq))
  nTrt = sum(samDataDF[, variable], na.rm = TRUE)
  nContr = nsamples(physeq) - nTrt
  if (min(nTrt, nContr) < 3.5 * nEval) {
    stop(
      "Not enough samples to make verification set two and a half times the size of the evaluation set!"
    )
  } else {
    
  }
  samNames = sample_names(physeq)
  samDataDF$set = rep("Not selected", nsamples(physeq))
  idEvalTrt = samNames[sample(which(samDataDF[, variable]), nEval)]
  idEvalControl = samNames[sample(which(!samDataDF[, variable]), nEval)]
  idVerifTrt = sample(samNames[!(samNames %in% idEvalTrt) &
                                 samDataDF[, variable]],
                      min(nTrt - nEval, maxVerif))
  idVerifControl = sample(samNames[!(samNames %in% idEvalControl) &
                                     !samDataDF[,
                                                variable]], min(min(nTrt, nContr) - nEval, maxVerif))
  samDataDF$set[samNames %in% idEvalTrt |
                  samNames %in% idEvalControl] = "eval"
  samDataDF$set[samNames %in% idVerifTrt |
                  samNames %in% idVerifControl] = "verif"
  samDataDF$group = ifelse(samDataDF[, variable], "grp1", "grp2")
  sample_data(physeq) = sample_data(samDataDF)
  # returnSeq = prune_samples(x=physeq, samDataDF$set!='Not selected')
  evalPhy = #simpleTrimGen(prune_samples(x = physeq, samDataDF$set == "eval"))
    verifPhy = #simpleTrimGen(prune_samples(x = physeq, samDataDF$set == "verif"))
    # Make sure only taxa present in both groups are retained
    taxaKeep = intersect(taxa_names(evalPhy), taxa_names(verifPhy))
  evalPhy = prune_taxa(x = evalPhy, taxaKeep)
  verifPhy = prune_taxa(x = verifPhy, taxaKeep)
  returnSeq = merge_phyloseq(evalPhy, verifPhy)
  rm(evalPhy, verifPhy)
  returnSeq
}
# 109 penicilin cases, 92 IBD cases

## Generate the Data
require(phyloseq, quietly = TRUE)
simListH1 <- lapply(simParamsH1[seq], function(iterRun) {
  set.seed(grep(iterRun, simParamsH1)) ## resets seed each iteration
  params <- strsplit(iterRun, delim)[[1]]
  names(params) <- simParamsLabelsH1
  
  ## write info about the current simulation on log-file sink(file =
  ## 'log0.txt', append = TRUE)
  cat(iterRun, "\t")
  # sink()
  
  # type of sample
  sampleTypeIter <- params["SampleType"]
  # The sample size to use for each group in this simulation
  nObs <- as.integer(params["nSamples"])
  # template and parameters
  template <- physeqList4Trim[[sampleTypeIter]]
  estPi <- piMoMs[[sampleTypeIter]]
  estTheta <- thetaMoMs[sampleTypeIter]
  estCov = covList[[sampleTypeIter]]
  estPhis = phiMLEs[[sampleTypeIter]]
  estRhos = rhoMLEs[[sampleTypeIter]]
  
  outLiers = OutLieList[[sampleTypeIter]]
  
  fC = as.numeric(params["EffectSize"])
  distrib = params["Distribution"]
  
  # Rarely a simulation has a weird value and fails.  Catch these with `try`,
  # and repeat the simulation call if error (it will be a new seed)
  tryAgain <- TRUE
  infLoopCount <- 1L
  maxLoops <- 15L
  
  while (tryAgain & infLoopCount <= maxLoops) {
    simResH1 <-
      microbioSim(
        postfix = iterRun,
        distrib = distrib,
        template = template,
        estPi = estPi,
        estTheta = estTheta,
        nObs = nObs,
        estRhos = estRhos,
        estPhis = estPhis,
        Covar = estCov,
        foldChange = fC,
        outLiers = outLiers,
        compensate = FALSE
      )
    
    ## Make sure there are at least 3 taxa per sample, even after trimming
    if (is.null(simResH1) | inherits(simResH1, "try-error")) {
      tryAgain <- TRUE
      infLoopCount <- infLoopCount + 1L
    } else {
      #simResH1 <- #simpleTrimGen(simResH1)
      if (fewTaxa(simResH1, 3L)) {
        tryAgain <- TRUE
        infLoopCount <- infLoopCount + 1L
      } else {
        tryAgain <- FALSE
      }  # END - ifelse: check not only error but *fewTaxa* success
    }  # END - ifelse: check successful simulation
  }  # END - while: protect against infinite loops and simulations failures
  
  if (infLoopCount > maxLoops) {
    warning("Consistent error found during simulation. Need to investigate cause.",
            immediate. = TRUE)
    cat(iterRun)
  } else {
    #simResH1 <- #simpleTrimGen(simResH1)
    
    ## only the _nOTUs_ most abundant OTUs are kept, unless the OTUs are already
    ## less
    if (ntaxa(simResH1) > nOTUs) {
      whichOTUs2Keep <- taxa_names(simResH1)[seq_len(nOTUs)]
      simResH1 <- prune_taxa(whichOTUs2Keep, simResH1)
    } else {
      
    }
  }  # END - ifelse: consistent error in current simulation
  
  ## log file writing
  cat("FINISHED\n")
  return(simResH1)
})  # END - parallelised simulations

#if (!file.exists("./brokenpromise/results/simulationListBetaBinomAH1.RData")) {
names(simListH1) <- simParamsH1[seq]
any(sapply(simListH1, class) != "phyloseq")
file_name <-  paste0("simulationListNegBin", letter, ".RData")
save(simListH1, simParamsLabelsH1, simParams, TPR, delim, file = file_name)
#   #} else {
#   #} else {
# #  }  # END - if: file with simulations under H1 already exists
# #}


## Matt's Code: Load/Annotate SimListH1
pr <- 0.1
BrokenPromiseData <- list()
for (index in 1:length(simListH1)) {
  physeq <- simListH1[[index]]
  sample_data <- sample_data(physeq) %>% as.data.frame
  sample_data$group <- c(rep("grp1", nrow(sample_data) / 2),
                         rep("grp2", nrow(sample_data) / 2)) %>% factor
  sample_data(physeq) <- sample_data(sample_data)
  otu_table(physeq) <- otu_table(physeq) %>%
    as("matrix") %>%
    #t %>%
    otu_table(taxa_are_rows = T)
  truede <- grep("-TP", rownames(otu_table(physeq)))
  params <- strsplit(names(simListH1)[index], "_")
  fx <- params[[1]][3]
  it <- params[[1]][2]
  sampleType <- params[[1]][1]
  distr <- params[[1]][5]
  iter <- paste0("|", it, "_", sampleType, "_", distr) ####
  m <- params[[1]][4]
  datalist <- list(physeq = physeq, truede = truede)
  BrokenPromiseData[[paste("BrokenPromiseH1", fx, m, pr, iter, sep = "_")]] <-
    datalist
}
distr <- ""
if (distribs == "negbinCorOut") {
  distr <- "NegBin"
} else if (distribs == "betabinCor") {
  distr <- "BetaBin"
}
                    
## the final simulated datasets to perform differential abundance tests on
save(BrokenPromiseData,
     file = paste0("BrokenPromiseDataJustSmallSamples_", distr, ".RData"))
