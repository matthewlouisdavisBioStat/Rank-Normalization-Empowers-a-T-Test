
      
     ## Make up an OTU Table with arbitrary values ##

otu_table <- matrix(rpois(100000,100)*rbinom(100000,1,0.2),nrow= 1000)
colnames(otu_table) <- sample(c("control sample","treatment sample"),
                              ncol(otu_table),replace= TRUE)
control_indices <- which(colnames(otu_table) == "control sample")
treatment_indices <- which(colnames(otu_table) == "treatment sample")
rownames(otu_table) <- paste0("OTU_",1:nrow(otu_table))



      ## Perform DAA using rank normalization with a t-test ## 

#1) Exclude OTU with only 0s
otu_table <- otu_table[rowSums(otu_table) > 0,] 

#2) Rank samples in ascending order, using average-ties
ranks <- apply(otu_table,2,rank) 

#3) For each OTU, perform t-test across ranks 
ttests <- apply(ranks,1,function(x)t.test(x[control_indices],
                                          x[treatment_indices]))
## Results for the first OTU
ttests$OTU_1








