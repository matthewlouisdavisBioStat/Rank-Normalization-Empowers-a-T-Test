## Matt's Function to Conduct Welch's Two-Sample T-Test Across All Rows at Once
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

