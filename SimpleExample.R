
  # Simple Demonstration of Rank Normalization Paired with a T-Test

## Make up an OTU Table with arbitrary values ## 
otu_table <- matrix(rpois(100000,100)*rbinom(100000,1,0.2),nrow= 1000)
colnames(otu_table) <- sample(c("control sample","treatment sample"),
                              ncol(otu_table),replace= TRUE)
control_indices <- which(colnames(otu_table) == "control sample")
treatment_indices <- which(colnames(otu_table) == "treatment sample")
rownames(otu_table) <- paste0("OTU_",1:nrow(otu_table))
head(otu_table[,1:5])

## Perform DAA using rank normalization with a t-test ## 

    #1) Exclude OTU with only 0s
otu_table <- otu_table[rowSums(otu_table) > 0,] 

    #2) Rank samples in ascending order, using average-ties
ranks <- apply(otu_table,2,rank) 

    #3) For each OTU, perform t-test across ranks 
results <- list() 
for(i in 1:nrow(ranks)){
  results[[paste0("OTU_",i)]] <- 
    t.test(ranks[i,control_indices],
           ranks[i,treatment_indices])
}

## Results for OTU 1 ##
results$OTU_1
