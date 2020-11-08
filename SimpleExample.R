
  # Simple demonstration of rank Normalization with a t-test using for loop

## Make up an OTU Table with Arbitrary Values
otu_table <- matrix(rnorm(100000),nrow= 1000)
colnames(otu_table) <- sample(c("control sample","treatment sample"),
                              ncol(otu_table),replace= TRUE)
control_indices <- which(colnames(otu_table) == "control sample")
treatment_indices <- which(colnames(otu_table) == "treatment sample")
rownames(otu_table) <- paste0("OTU_",1:nrow(otu_table))
head(otu_table[,1:5])

## Perform DAA using Rank Normalization with a t-test
otu_table <- otu_table[rowSums(otu_table) > 0,] #1) Exclude OTU with only 0s
ranks <- apply(otu_table,2,rank) #2) Rank samples in ascending order
results <- list() #3) Perform t-test across ranks with respect to the two groups
for(i in 1:nrow(ranks)){
  results[[paste0("OTU_",i)]] <- 
    t.test(ranks[i,control_indices],
           ranks[i,treatment_indices])
}

## Results for OTU 1
results$OTU_1


################################################################################

    ## Increase the speed of a t-test with vectorization

## Matt's function to conduct Welch's Two-Sample t-test across all rows at once
fastT <- function(otu, ## otu_table, matrix with rows as taxa columns as samples
                  indg1, ## treatment indices
                  indg2){ ## control indices
  l1 <- length(indg1)
  l2 <- length(indg2)
  S1 <- (rowSums((otu[,indg1] - rowMeans(otu[,indg1]))^2) / (l1 - 1))
  S2 <- (rowSums((otu[,indg2] - rowMeans(otu[,indg2]))^2) / (l2 - 1))
  sigmas <- sqrt ( S1/l1 + S2/l2 )
  df <- sigmas^4  / ( ((l1^2 * (l1-1))^-1 *  S1^2) + ((l2^2 * (l2-1))^-1 *  S2^2) )
  ## 
  obsT <- (rowMeans(otu[,indg1]) - rowMeans(otu[,indg2])) / sigmas
  ## 
  pvals <- 2*(1-pt(abs(obsT), df = df))
  ifelse(is.nan(pvals), 1, pvals)
}

## returns raw p-values
results <- fastT(ranks,treatment_indices,control_indices)
head(results)
