
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
minReads <- 1
minPrev <- 0.05
prevalence <- rowMeans(otu_table >= minReads)
indOTUs2Keep <- (prevalence >= minPrev)
nonzeroind <- indOTUs2Keep

  ## Rank Normalization with a T-Test
otu_table <- otu_table[rowSums(otu_table) > 0,] # 
ranks <- apply(otu_table, 2, rank)
l1 <- length(indg1)
l2 <- length(indg2)
S1 <- (rowSums((ranks[,indg1] - rowMeans(ranks[,indg1]))^2) / (l1 - 1))
S2 <- (rowSums((ranks[,indg2] - rowMeans(ranks[,indg2]))^2) / (l2 - 1))
sigmas <- sqrt ( S1/l1 + S2/l2 )
df <- sigmas^4  / ( ((l1^2 * (l1-1))^-1 *  S1^2) + ((l2^2 * (l2-1))^-1 *  S2^2) )
obsT <- (rowMeans(ranks[,indg1]) - rowMeans(ranks[,indg2])) / sigmas
pvals <- 2*(1-pt(abs(obsT), df = df))
pvals <- ifelse(is.nan(pvals), 1, pvals)
names(pvals) <- rownames(ranks)

# padj <- rep(1, length(pvals))
pvals <- pvals[names(nonzeroind)]
padj <- p.adjust(pvals, method = "BH")
names(padj) <- names(pvals)

mat <- cbind(c(pvals), c(padj))
colnames(mat) <- c("rawP", "adjP")
rownames(mat) <- names(pvals)

test <- "t-test"
normFacts <- "rank-simple"
reject <-  names(padj)[which(padj < .05)]
degenes <- degenes

cor_rej <- sum(reject %in% degenes)
err_rej <- sum(reject %!in% degenes)
Specificity <- 1-err_rej / ( nrow(otu_table(physeq)) -length(degenes) )
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
tTest_rank  <- list("res" = res,
                    "pvals" = rawP)
Sys.time() - start

##############################################################
