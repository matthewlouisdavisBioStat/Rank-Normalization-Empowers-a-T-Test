
#####################################################
  ## Load Data
start <- Sys.time()
otu_table <- as(otu_table(physeq), "matrix")
ncol <- ncol(otu_table)
nrow <- nrow(otu_table)
sample_data <- as(sample_data(physeq), "data.frame")
indg1 <- which(sample_data$group == "grp1")
indg2 <- which(sample_data$group == "grp2")

  ## Remove rows containing only 0s, but make note of the 5% prevalant OTU
  ## Recall: Other methodologies prefer 5% prevalance 

  ## Rank Normalization with a T-Test
otu_table <- otu_table[rowSums(otu_table) > 0,] # 
minReads <- 1
minPrev <- 0.05
prevalence <- rowMeans(otu_table >= minReads)
indOTUs2Keep <- (prevalence > minPrev)
nonzeroind <- indOTUs2Keep
ranks <- apply(otu_table, 2, rank)
#ranks <- ranks[names(nonzeroind),]


  ## function worked faster if I broke up into 10 little small parts
pvals <- rowMeans(sapply(1:10,function(x){
  #set.seed(x)
  permTest(ranks, 
                  N = 1000,
                  indg1 = indg1,
                  indg2 = indg2,
                  makeCluster = FALSE,
                  adj = "")})) 

#,
pvals <- pvals[names(nonzeroind)]
#padj <- pvals
padj <- p.adjust(pvals,"BH")
# #head(pvals)
# names(pvals) <- rownames(ranks)
# 
# # padj <- rep(1, length(pvals))
# pvals <- pvals[names(nonzeroind)]
# padj <- p.adjust(pvals, method = "BH")
# names(padj) <- names(pvals)
# 
mat <- cbind(c(pvals), c(padj))
colnames(mat) <- c("rawP", "adjP")
rownames(mat) <- names(pvals)

test <- "perm-test"
normFacts <- "rank"
reject <-  names(padj)[which(padj < .05)]

degenes <- degenes

cor_rej <- sum(reject %in% degenes)
err_rej <- sum(reject %!in% degenes)
Specificity <- err_rej / ( nrow(otu_table(physeq)) -length(degenes) )
Sensitivity <- cor_rej / length(degenes)
FDR <- err_rej / (cor_rej + err_rej)

res <- matrix(c( test,
                 normFacts,
                 Sensitivity,
                 Specificity,
                 FDR ), nrow = 1)
colnames(res) <- c( "test",
                    "norm",
                    "sens",
                    "spec",
                    "fdr" )
res[,"fdr"] <- ifelse(is.nan(res[,"fdr"]), 0, res[,"fdr"])
res[,"fdr"] <- ifelse(res[,"fdr"] == "NaN", 0, res[,"fdr"])
rawP <- pvals
permTest_rank  <- list("res" = res,
                    "pvals" = rawP)
Sys.time() - start

##############################################################
