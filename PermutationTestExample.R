
# Simple Demonstration of Rank Normalization With A Permutation T-Test

################################################################################

## Nonparametric Permutation T-Test ##
permTest <- function(    otu_table,    # rows as taxa, columns as samples
                         ncores = 8, # number of cores used
                         N = 1000,   # number of permutations
                         indg1,      # column indices of group 1
                         indg2,      # column indices of group 2
                         alpha = 0.05, # significance level 
                         adj = "",   # BH? FDR? what kind of adjustment, if any
                         makeCluster = T, # if we should make a cluster within the function (unless we make one outside)
                         juststat = FALSE){ # return only the permuted t-statistics
  ## Observed T-Statistic
  l1 <- length(indg1)
  l2 <- length(indg2)
  S1 <- (rowSums((otu_table[,indg1] - rowMeans(otu_table[,indg1]))^2) / 
           (l1 - 1))
  S2 <- (rowSums((otu_table[,indg2] - rowMeans(otu_table[,indg2]))^2) / 
           (l2 - 1))
  sigmas <- sqrt ( S1/l1 + S2/l2 )
  TObs <- (rowMeans(otu_table[,indg1]) - rowMeans(otu_table[,indg2])) / 
            sigmas
  TObs <- ifelse(is.nan(TObs), 0, TObs)
  
  ## Run in Parallel
  ncol <- ncol(otu_table)
  if(makeCluster){
    cl <- makeCluster(ncores)
    doSNOW::registerDoSNOW(cl)
  }
  TPerm <- foreach(i = 1:N, .combine = "cbind") %dopar% {
    ## Permute Column Labels
    otu_tablePerm <- otu_table[,sample(1:ncol, ncol, replace = FALSE)]
    ## Recalc Statistic
    S1 <- (rowSums((otu_tablePerm[,indg1] - 
                      rowMeans(otu_tablePerm[,indg1]))^2) / (l1 - 1))
    S2 <- (rowSums((otu_tablePerm[,indg2] - 
                      rowMeans(otu_tablePerm[,indg2]))^2) / (l2 - 1))
    sigmas <- sqrt ( S1/l1 + S2/l2 )
    TP <- (rowMeans(otu_tablePerm[,indg1]) - 
             rowMeans(otu_tablePerm[,indg2])) / sigmas
    ifelse(is.nan(TP), 0, TP)
  }
  if(makeCluster){
    stopCluster(cl)
  }
  if(juststat == TRUE){
    return(TPerm)
  }
  
  if(adj == "BH"){
    ## BH correction for pvalues 
    return(rownames(TPerm)[which(p.adjust(rowMeans(abs(TPerm) > abs(TObs)),"BH") <= alpha)])
  } else if(adj == "fdr"){
    ## FDR stepdown correction procedure for controlling pvalues
    cs <- seq(0.00001,max(abs(TObs))-.00001, length.out = 10000)
    fdrs <- rep(NA, 10000)
    names(fdrs) <- cs
    if(makeCluster){
      cl <- makeCluster(4)
      doSNOW::registerDoSNOW(cl)
    }
    fdrs <- foreach(i = 1:length(cs),.combine = "c") %dopar% {
      c <- cs[i]
      R <- sum(abs(TObs) >= c)
      W <- mean(colSums(TPerm >= c))
      W/R
    }
    if(makeCluster){
      stopCluster(cl)
    }
    diff <- abs(fdrs - alpha)
    names(diff) <- cs
    C <- as.numeric(names(diff)[which(diff == min(diff, na.rm = TRUE))])
    ## rejections
    return(names(TObs)[which(abs(TObs) > abs(C))])
  } else{
    ## Return raw pvalues otherwise
    return(rowMeans(abs(TPerm) > abs(TObs)))
  }
  
}
################################################################################


## Create an OTU Table with Arbitrary Values
otu_table <- matrix(rpois(100000,100),nrow= 1000)
colnames(otu_table) <- sample(c("control sample","treatment sample"),
                              ncol(otu_table),replace= TRUE)
control_indices <- which(colnames(otu_table) == "control sample")
treatment_indices <- which(colnames(otu_table) == "treatment sample")
rownames(otu_table) <- paste0("OTU_",1:nrow(otu_table))
head(otu_table[,1:5])

## Rank normalization
otu_table <- otu_table[rowSums(otu_table) > 0,]
ranks <- apply(otu_table,2,rank)

## Run Permutation t-test non-parametric alternative to Welchs t-test
## Returns pvals of each OTU by default
library(foreach)
library(doSNOW)
library(parallel)

cl <- makeCluster(8)
registerDoSNOW(cl)
pvals <- permTest(ranks, 
                N = 500,
                indg1 = control_indices,
                indg2 = treatment_indices) #,
stopCluster(cl)
pvals
