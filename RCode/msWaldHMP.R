

    ## Matt: I copied and pasted this from Mr. Mattiello's github
    ## The Github download wasn't working

#msWaldFunctions
################################################################################
### Title: internals.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 04/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 04/mar/2015 - 16:21:48:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################

###### my version of Weir computation (Method of Moments)
weirMoM4Wald <- function (data, se = FALSE)
{
  if (missing(data)) 
    stop("data missing.")
  K <- ncol(data)
  J <- nrow(data)
  
  totalsample <- sum(data, na.rm = TRUE)
  MoM <- colSums(data, na.rm = TRUE)/totalsample
  Sn <- rowSums(data, na.rm = TRUE)
  auxData <- data / Sn
  
  MSP <- sum(rowSums(
    (auxData - matrix(rep(MoM, J), nrow = J, ncol = K, byrow = TRUE))^2,
    na.rm = TRUE) * Sn) / (J - 1)
  MSG <- sum(rowSums(auxData * (1 - auxData), na.rm = TRUE) * Sn) / (totalsample - J)
  
  nc <- (sum(Sn) - sum(Sn^2)/sum(Sn)) / (J - 1)
  MoM.wh <- (MSP - MSG)/(MSP + (nc - 1) * MSG)
  
  if (se)
  {
    std.er <- sqrt(2 * (1 - MoM.wh)^2 / (J - 1) * 
                     ((1 + (nc - 1) * MoM.wh)/nc)^2)
    return(list(theta = MoM.wh, se = std.er))
  } else
  {
    return(MoM.wh)
  }
}# END - myWeirMom


### pi vector estimation with Method of Moments
piMoM4Wald <- function(data)
{
  totalReads <- sum(data, na.rm = TRUE)
  piMom <- colSums(data, na.rm = TRUE)/totalReads
  zeroInds <- abs(piMom) < .Machine$double.eps
  r <- sum(zeroInds)
  rr <- length(piMom) - r
  piMom[!zeroInds] <- piMom[which(piMom != 0)] - r/(rr * 2 * (totalReads + 1))
  piMom[zeroInds] <- 1/(2 * (totalReads + 1))
  
  return(piMom)
}

################################################################################
### Title: msWald.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 04/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 04/mar/2015 - 16:27:52:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Wrapper for msWaldStat
#' 
#' @title Wrapper for msWaldStat
#' 
#' @param nReads 
#'   either (i) a \code{list} or (ii) a \code{list} of \code{list}s; more in detail:
#'   \itemize{
#'     \item{}{
#'     in case (i) \code{nReads} has to be a list of number of reads (library size) for
#'     each group, {\em i.e.} the i-th element of the list contains \eqn{n_i} elements;
#'     }
#'     \item{}{
#'     in case (ii) \code{nReads} must be a \code{list} where each element contains a 
#'     \code{list} with the same structure of the previous point, {\em i.e.} each element
#'     element is used to perform a separate analysis on the same data by subsetting the
#'     samples. It useful in case of simulations when several sample sizes are to be 
#'     tried.
#'     }
#'   }
#'   See also the {\em details} section.
#' @param alphaDM 
#'   \eqn{\alpha}{alpha} parameter for the Dirichlet-Multinomial distribution.
#'   Either a \code{matrix} of dimensions \eqn{K \times G}{K x G}, K = number of OTUs, 
#'   G = number of groups, containing the parameter vectors for the Dirichlet 
#'   distribution or a \code{list} where each element contains a \eqn{K \times G}{K x G}
#'   \code{matrix} in case the simulations need to be stratified, {\em e.g.} by 
#'   enterotype. In the latter case simulations are performed separately 
#'   stratum-by-stratum; individual and global rejection rates (power) are given as 
#'   output (see also details).
#' does NOT need to be stratified (subsets of samples, 
#'     {\em e.g.} in case of enterotypes)
#' @param thetaDM 
#'   \eqn{\theta}{theta} overdispersion parameter for the Dirichlet-Multinomial 
#'   distribution. Either (i) a \code{numeric} vector of length equal 
#'   to the number of groups under test or (ii) a \code{matrix} where each column 
#'   corresponds to a different stratum. It is recycled when possible (same values for 
#'   each stratum).
#' @param wmnTest 
#'   \code{logical} default to \code{FALSE}, if \code{TRUE} performs the two-samples 
#'   Wilcoxon-Mann-Whitney test one OTU at-a-time and then corrects for multiplicity.
#'   {\bfseries Note:} it works only with case (i), thus without subsets and 
#'   without strata.
#' @param adjMethod
#'   \code{character} method to adjust for multiplicity in case \code{wmnTest} is 
#'   \code{TRUE}
#' @return 
#'     what the function returns
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
### nReads is a list of number of reads, i-th position of the list contains 
### n_i elements: number of reads of that observation/sample (library size)
msWald <- function(nReads, alphaDM, thetaDM, wmnTest = FALSE, adjMethod = "fdr")
{
  ### check class and names of *nReads*
  if (!is.list(el(nReads)))
  {
    nReads <- list("totSet" = nReads)
  } else
  {
    if (is.null(names(nReads)))
      names(nReads) <- paste0("subset", seq_along(nReads))
  }
  
  ### check class and names of *alphaList*
  if (!is.list(alphaDM))
  {
    alphaDM <- list(alphaDM)
  } else
  {
    if (is.null(names(alphaDM)))
      names(alphaDM) <- paste0("stratum", seq_along(alphaDM))
  }
  
  nGroups <- length(el(nReads))
  nStrata <- length(alphaDM)
  nSubsets <- length(nReads)
  
  ### checks of thetaVec
  thetaDM <- as.matrix(thetaDM)
  if (ncol(thetaDM) < nStrata)
  {
    thetaDM <- matrix(thetaDM, nrow = length(thetaDM), ncol = nStrata)
    #    message("  theta parameters has been recycled")
  } else {}
  rownames(thetaDM) <- paste0("group", seq_len(nGroups))
  colnames(thetaDM) <- names(alphaDM)
  
  
  sampleSizes <- sapply(nReads, FUN = function(x) sapply(x, FUN = length))
  rownames(sampleSizes) <- paste0("group", seq_len(nGroups))
  colnames(sampleSizes) <- names(nReads)
  maxSampleSize <- apply(sampleSizes, MARGIN = 1, max)
  
  dmDataList <- vector("list", nGroups)
  names(dmDataList) <- paste0("group", seq_len(nGroups))
  
  
  piMomAux <- array(NA, dim = c(dim(el(alphaDM)), nSubsets), 
                    dimnames = list(
                      rownames(el(alphaDM)), colnames(el(alphaDM)), 
                      names(nReads))
  )
  
  thetaMomAux <- matrix(NA, nrow = nGroups, ncol = nSubsets, 
                        dimnames = list(
                          paste0("group", seq_len(nGroups)), 
                          names(nReads)
                        ))
  
  res <- matrix(NA, nrow = nStrata, ncol = nSubsets, 
                dimnames = list(names(alphaDM), names(nReads)))
  
  
  aux <- nReads[apply(sampleSizes, 1, which.max)]
  nReadsMax <- lapply(seq_along(aux), 
                      FUN = function(i, dat) dat[[c(i, i)]], dat = aux)
  
  
  ### main loop over strata (enterotypes), usually equal to 1 or 3
  for (stRun in seq_len(nStrata))
  {
    ## loop over groups (can be more than 2)
    for (grRun in seq_len(nGroups))
    {
      dmDataList[[grRun]] <- t(sapply(seq_len(maxSampleSize[grRun]), 
                                      FUN = function(i, nRds, alMat)
                                      {
                                        rmultinom(n = 1, size = nRds[i], prob = rDirichlet(n = 1, alpha = alMat))
                                      }, 
                                      nRds = nReadsMax[[grRun]], 
                                      alMat = alphaDM[[stRun]][, grRun] * 
                                        (1 - thetaDM[grRun, stRun])/thetaDM[grRun, stRun]
      ))
      
      ## populate *theta* and *pi* for the Wald statistic
      for (subRun in seq_len(nSubsets))
      {
        subsetData <- dmDataList[[grRun]][seq_len(sampleSizes[grRun, subRun]), ]
        thetaMomAux[grRun, subRun] <- msWaldHMP:::weirMoM4Wald(subsetData)
        piMomAux[, grRun, subRun] <- msWaldHMP:::piMoM4Wald(subsetData)
      }# END - for: subRun
    }# END - for: grRun
    
    ## compute Wald statistics for each subset in the stratum
    for (subRun in seq_len(nSubsets))
    {
      res[stRun, subRun] <- msWaldHMP:::msWaldStat(
        nReads[[subRun]], piMat = piMomAux[, , subRun], thetaVec = thetaMomAux[, subRun])
    }# END - for: subRun
  }# END - for: stRun 
  
  
  
  if (length(el(nReads)) == 2L && wmnTest)
  {
    nReads <- nReads[[1L]]
    alphaDM <- alphaDM[[1L]]
    
    pVals <- sapply(seq_len(nrow(alphaDM)), FUN = function(i, dat1, dat2)
    {
      wilcox.test(x = dat1[, i], y = dat2[, i], exact = FALSE)$p.value 
    }, 
    dat1 = dmDataList[[1L]], dat2 = dmDataList[[2L]]
    )
    
    pValsAdj <- p.adjust(pVals, method = adjMethod)
    res <- list(out, cbind("pVals" = pVals, "pValsAdj" = pValsAdj))
  } else {}
  
  res
  
}# END - msWald


################################################################################
### Title: msWaldMC.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 04/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 04/mar/2015 - 16:31:08:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Wrapper function for MonteCarlo simulations
#' 
#' @title MonteCarlo Simulations with Multi-Sample Wald Statistic
#' 
#' @param MC 
#'   number of MonteCarlo iterations
#' @param nReads 
#'   either (i) a \code{list} or (ii) a \code{list} of \code{list}s; more in detail:
#'   \itemize{
#'     \item{}{
#'     in case (i) \code{nReads} has to be a list of number of reads (library size) for
#'     each group, {\em i.e.} the i-th element of the list contains \eqn{n_i} elements;
#'     }
#'     \item{}{
#'     in case (ii) \code{nReads} must be a \code{list} where each element contains a 
#'     \code{list} with the same structure of the previous point, {\em i.e.} each element
#'     element is used to perform a separate analysis on the same data by subsetting the
#'     samples. It useful in case of simulations when several sample sizes are to be 
#'     tried.
#'     }
#'   }
#'   See also the {\em details} section.
#' @param alphaDM 
#'   \eqn{\alpha}{alpha} parameter for the Dirichlet-Multinomial distribution.
#'   Either a \code{matrix} of dimensions \eqn{K \times G}{K x G}, K = number of OTUs, 
#'   G = number of groups, containing the parameter vectors for the Dirichlet 
#'   distribution or a \code{list} where each element contains a \eqn{K \times G}{K x G}
#'   \code{matrix} in case the simulations need to be stratified, {\em e.g.} by 
#'   enterotype. In the latter case simulations are performed separately 
#'   stratum-by-stratum; individual and global rejection rates (power) are given as 
#'   output (see also details).
#' does NOT need to be stratified (subsets of samples, 
#'     {\em e.g.} in case of enterotypes)
#' @param thetaDM 
#'   \eqn{\theta}{theta} overdispersion parameter for the Dirichlet-Multinomial 
#'   distribution. Either (i) a \code{numeric} vector of length equal 
#'   to the number of groups under test or (ii) a \code{matrix} where each column 
#'   corresponds to a different stratum. It is recycled when possible (same values for 
#'   each stratum).
#' @param sigLev 
#'   significance level (alpha, or type-I error)
#' @param avgRej 
#'   if FALSE it returns the number of rejections instead of the proportion 
#'   (among MC iterations)
#' @param ...
#'   see \code{\link{msWald}} for details
#' @return 
#'     either number of rejections or the mean
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
msWaldMC <- function(MC = 100, nReads, alphaDM, thetaDM, sigLev = .05, avgRej = TRUE,
                     wmnTest = FALSE, adjMethod = "fdr")
{
  ### check class and names of *nReads*
  if (!is.list(el(nReads)))
  {
    nReads <- list("totSet" = nReads)
  } else
  {
    if (is.null(names(nReads)))
      names(nReads) <- paste0("subset", seq_along(nReads))
  }
  
  ### check class and names of *alphaList*
  if (!is.list(alphaDM))
  {
    alphaDM <- list(alphaDM)
  } else
  {
    if (is.null(names(alphaDM)))
      names(alphaDM) <- paste0("stratum", seq_along(alphaDM))
  }
  
  nGroups <- length(el(nReads))
  nStrata <- length(alphaDM)
  nSubsets <- length(nReads)
  nOtus <- sapply(alphaDM, FUN = nrow)
  
  ### quantiles of the reference distribution
  ## Degrees of Freedom
  DoFs <- (nGroups - 1) * (nOtus - 1)
  qAlpha <- qchisq(p = 1 - sigLev/nStrata, df = DoFs, ncp = 0, lower.tail = TRUE)
  qAlphaGlob <- qchisq(p = 1 - sigLev, df = sum(DoFs), ncp = 0, lower.tail = TRUE)
  
  ### MonteCarlo simulations
  tmp <- lapply(seq_len(MC), 
                FUN = function(x, wmn, adj)
                {
                  msWald(nReads, alphaDM, thetaDM, wmnTest = wmn, adjMethod = adj)
                }, wmn = wmnTest, adj = adjMethod)
  
  res <- array(unlist(tmp), dim = c(nStrata, nSubsets, MC), 
               dimnames = list(
                 rownames(el(tmp)), colnames(el(tmp)), paste0("MC", seq_len(MC))
               ))
  
  
  ### TODO: use *pchisq* plus multiplicity correction instead of just checking rejection
  #  resAdj <- pchisq(res, df = DoFs[1L], lower.tail = FALSE)
  #  resAdj <- apply(resAdj, MARGIN = 2L:3, fun = p.adjust, method = adjMethod)
  #  rejRes <- rowSums(resAdj <= sigLev, na.rm = TRUE, dims = 2)
  
  
  ### number of rejections, note that the following command is equivalent to:
  ### rejRes <- apply(res > qAlpha, MARGIN = 1L:2, FUN = sum) 
  rejRes <- rowSums(res > qAlpha, na.rm = TRUE, dims = 2)
  
  if (avgRej)
  {
    rejRes <- rejRes / MC
  } else {}# END - ifelse: avgRej
  
  
  ### if multilple strata, sum up together
  if (nStrata > 1L)
  {
    ##   TODO: sums wald stat and compare with different quatile of chisquare
    #    rejResGlob <- rowSums(apply(res > qAlpha, MARGIN = c(2L, 3), 
    #        FUN = function(x) !all(!x)))
    rejResGlob <- rowSums(colSums(res, na.rm = TRUE) > qAlphaGlob)
    
    if (avgRej)
    {
      rejResGlob <- rejResGlob / MC
    }
    
    ## bind together the two
    rejRes <- rbind(rejRes, "Global" = rejResGlob)
  } else {}
  
  
  
  
  
  
  if (!is.null(dim(rejRes)) && wmnTest)
  {
    waldRes <- unlist(res[2 * seq_len(MC) - 1L]) > qAlpha
    
    wmnRes <- sapply(res[2 * seq_len(MC)], function(x, al)
    {
      colSums(x <= al, na.rm = TRUE) > 0
    }, al = sigLev)
    
    if(avgRej)
    {
      c("wald" = mean(waldRes, na.rm = TRUE), "WMN" = rowMeans(wmnRes, na.rm = TRUE))
    } else
    {
      c("wald" = sum(waldRes, na.rm = TRUE), "WMN" = rowSums(wmnRes, na.rm = TRUE))
    }# END - ifelse: avgRej
  }# END - ifelse: wmnTest
  
  
  rejRes
  
}

#set.seed(12345)
#ha <- msWaldMC(MC = 500, nReads = nReads, alphaDM = alphaDM, thetaDM = thetaDM)
#set.seed(12345)
#ha2 <- msWaldMC(MC = 500, nReads = nReads[[3L]], alphaDM = alphaDM, thetaDM = thetaDM)
#set.seed(12345)
#ha3 <- msWaldMC(MC = 500, nReads = nReads, alphaDM = alphaDM$ent1, thetaDM = thetaDM[, 1])
#
#ha
#ha2
#ha3


################################################################################
### Title: msWaldStat.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 04/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 04/mar/2015 - 16:22:38:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Multi-sample Wald-type test statistics, basically Xmcupo.sevsample()
#' 
#' @title Multi-Sample Wald Statistic
#' 
#' @param nReads 
#'   list of number of reads (library size) for each group, 
#'   i-th element of the list contains \eqn{n_i} elements
#' @param piMat 
#'   matrix K x G, K = number of OTUs, G = number of groups, contains probability 
#'   vectors/R.A.D.s
#' @param thetaVec 
#'   overdispersion parameter vector, one for each group
#' @return 
#'     \code{numeric} Wald-type test statistic
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
msWaldStat <- function(nReads, piMat, thetaVec)
{
  nGroups <- length(nReads)
  
  N_1Cj <- rep.int(NA, nGroups)
  for (grRun in seq_len(nGroups))
  {
    N_1Cj[grRun] <- (thetaVec[grRun] * (sum(nReads[[grRun]]^2) - sum(nReads[[grRun]])) + 
                       sum(nReads[[grRun]])) / sum(nReads[[grRun]])^2
  }
  
  N_1CjInv <- 1 / N_1Cj
  
  pi0Est <- drop((piMat %*% N_1CjInv) / sum(N_1CjInv))
  
  Xmcupo <- sum(N_1CjInv * colSums(
    (piMat - matrix(pi0Est, nrow = length(pi0Est), ncol = nGroups))^2 / drop(pi0Est)
  ))
  
  Xmcupo
}# END - msWaldStat

################################################################################
### Title: rDirichlet.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 04/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 04/mar/2015 - 16:20:08:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Draws from Dirichlet distribution
#' 
#' @title Random Dirichlet Draws
#' 
#' @param n
#'     number of draws
#' @param alpha
#'     parameter vector
#' @return 
#'     what the function returns
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 

rDirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}


################################################################################
### Title: rDirichlet.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 04/mar/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 04/mar/2015 - 16:20:08:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Draws from Dirichlet distribution
#' 
#' @title Random Dirichlet Draws
#' 
#' @param n
#'     number of draws
#' @param alpha
#'     parameter vector
#' @return 
#'     what the function returns
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 

rDirichlet <- function (n, alpha) 
{
  l <- length(alpha)
  x <- matrix(rgamma(l * n, alpha), ncol = l, byrow = TRUE)
  sm <- x %*% rep(1, l)
  x/as.vector(sm)
}


################################################################################
### Title: wrapperWMW.R                                              ###
###                                                                          ###
### Project: msWaldHMP                                                ###
###                                                                          ###
### Version: 0.1 - 08/set/2015                                    ###
###                                                                          ###
### Description: short description                                           ###
###                                                                          ###
### Authors: Federico Mattiello <Federico.Mattiello@UGent.be>                ###
###                                                                          ###
### Maintainer: Federico Mattiello <Federico.Mattiello@UGent.be>             ###
###                                                                          ###
### Versions:                                                                ###
###     > 0.1 - 08/set/2015 - 14:13:31:                                      ###
###         creation                                                         ###
###                                                                          ###
################################################################################
###
#' Wrapper function for Wilcoxon-Mann-Whitney test
#' 
#' Utility function performing WMW test on each OTU individually, and then correcting 
#' for multiplicity with one of the standard corrections, \code{"fdr"} is the default. 
#'
#' @title Wrapper function for Wilcoxon-Mann-Whitney test
#' 
#' @param x
#'     \code{matrix} containing the first sample
#' @param y
#'     \code{matrix} containing the second sample
#' @param ...
#'     other parameters for \link{\code{wilcox.test}}
#' @return 
#'     \code{list} with two elements of length \code{ncol(x)}, 
#'     \code{"statistic"} containing the test-statistics, and 
#'     \code{"p.value"} containing the adjusted p.values vector.
#' @author Federico Mattiello <Federico.Mattiello@@UGent.be>
#' @export
#' 
wrapperWMW <- function(x, y, adjMethod = "BH", ...)
{
  res <- lapply(seq_len(NCOL(x)), 
                FUN = function(i)
                {
                  suppressWarnings(stats:::wilcox.test.default(
                    x = x[, i], y = y[, i], ...)[c("statistic", "p.value")])
                  
                })
  rawPvals <- sapply(res, elNamed, name = "p.value")
  
  list(
    "statistic" = sapply(res, elNamed, name = "statistic"),
    "p.value" = p.adjust(rawPvals, method = adjMethod)
  )
}# END: function - wrapperWMW.R








