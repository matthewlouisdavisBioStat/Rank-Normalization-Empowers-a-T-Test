
# Simple Demonstration of Rank Normalization with a T-Test

## Make up an OTU Table with Arbitrary Values
otu_table <- matrix(rpois(100000,100),nrow= 1000)
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
