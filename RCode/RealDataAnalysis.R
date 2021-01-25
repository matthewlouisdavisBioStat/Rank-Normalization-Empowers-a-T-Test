
## Real Data Analysis ##

## Note: I'm using R 3.6.3 Right Now ##

## Choose Your Data, Zeller or RISK

## Download xlsx files directly from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4299606/,
#     "msb0010-0766-sd4.xlsx"
#     "msb0010-0766-sd2.xlsx"
## and save in datawd

## Load MicrobeDS curated data
dataset <- c("RISK"); library(MicrobeDS)
#dataset <- c("Zeller"); datawd <- "C:/Users/Matthew/Documents/Courses/Kai/Final Results/RCode/"

## Load Libraries
library(dplyr)
library(readxl)
library(metagenomeSeq)
library(edgeR)
library(limma)
library(DESeq2)
library(ALDEx2)
library(phyloseq)
################################################################################

## Helper Functions


'%!in%' <- function(x,y)!(x%in%y)

set_colnames <- function(x,colnames){
  colnames(x) <- colnames
  x
}

## Performs T-Test Simultaneously Across OTU Table
fastT <- function(otu_table, indg1, indg2, sign = F,stat = "no"){
  l1 <- length(indg1)
  l2 <- length(indg2)
  S1 <- (rowSums((otu_table[,indg1] - rowMeans(otu_table[,indg1]))^2) / (l1 - 1))
  S2 <- (rowSums((otu_table[,indg2] - rowMeans(otu_table[,indg2]))^2) / (l2 - 1))
  sigmas <- sqrt ( S1/l1 + S2/l2 )
  df <- sigmas^4  / ( ((l1^2 * (l1-1))^-1 *  S1^2) + ((l2^2 * (l2-1))^-1 *  S2^2) )
  obsT <- (rowMeans(otu_table[,indg1]) - rowMeans(otu_table[,indg2])) / sigmas
  if(stat == "yes"){
    return(obsT)
  } else{
    pvals <- 2*(1-pt(abs(obsT), df = df))
    if(sign){
      pvals <- pvals*sign(obsT)
    }
    ifelse(is.nan(pvals), 1, pvals)
  }
}

################################################################################

## Load Data

## Import and Process Data of Zeller
if(dataset == "Zeller"){
  otu_table <- readxl::read_xlsx(paste0(datawd,"msb0010-0766-sd4.xlsx"))
  otu_table <- as.data.frame(otu_table[2:nrow(otu_table),])
  rownames(otu_table) <- otu_table[,1]
  otu_table <- otu_table[,2:ncol(otu_table)]
  readData <- readxl::read_xlsx(paste0(datawd, "msb0010-0766-sd2.xlsx"), sheet = 3)
  rawReads <- readData$`Raw Reads`
  names(rawReads) <- readData$'Sample ID'
  rawReads <- rawReads[which(names(rawReads) %in% colnames(otu_table))]
  rawReads <- as.numeric(rawReads)
  for(j in 1:ncol(otu_table)){
    otu_table[,j] <- otu_table[,j]*rawReads[j]
  }
  sample_data <- readxl::read_xlsx(paste0(datawd,"msb0010-0766-sd2.xlsx"),
                                   sheet = 2) %>% as.data.frame %>%
    subset(Diagnosis %in% c("Normal","Cancer")) %>%
    as.data.frame
  
  rownames(sample_data) <- sample_data$`Sample ID`
  sample_data <- sample_data[colnames(otu_table),]
  otu_table <- as(otu_table,"matrix")
  otu_table <- otu_table[,colnames(otu_table) %in% rownames(sample_data)]
  sample_data <- sample_data[colnames(otu_table),]
  indg1 <- which(sample_data$Group == "CRC")
  indg2 <- which(sample_data$Group == "Control")
  otu_table <- otu_table[rowSums(otu_table)>0,]
  design <- model.matrix(as.formula("~ Group"), data = as.data.frame(as(sample_data,"matrix")))
}

## Test if Demographic Variables May Be Controlled For

# fit <- glm(factor(sample_data$Group) ~ sample_data$Gender, 
# family = "binomial") # Alt Model
# coef(summary(fit))
# 
# fit <- glm(factor(sample_data$Group) ~ 1, 
# family = "binomial") # Null Model
# coef(summary(fit))


## Microbe RISK CCFA dataset
if(dataset == "RISK"){
  physeq <- MicrobeDS::RISK_CCFA
  tax <- tax_table(physeq) %>% as("matrix") %>% as.data.frame
  sample_data <- sample_data(physeq) %>% as.data.frame %>%
    subset(antibiotics == "false" & ## no antibiotics
             collection == "RISK" & ## isolates only the RISK patients, not other studies
             diagnosis %in% c("CD","no") & ## two disease statuses of interest, excludes UC etc.
             biopsy_location %in% c("Terminal ileum") & ## the most common biopsy location
             diseasesubtype %in% c("iCD","no") & ## most common type of CD
             sample_type == "biopsy" & ## Only 
             b_cat == "B1" & ## Removed B2,B3
             
             ## Remove subjects with symptoms of other conditions
             !(inflammationstatus == "inflamed" & diagnosis == "no") &
             !(ileal_invovlement == "true" & diagnosis == "no") & 
             !(gastric_involvement == "true" & diagnosis == "no") & 
             !(disease_extent != "missing" & diagnosis == "no") & 
             !(perianal != "false" & diagnosis == "no") &
             
             ## Isolate females 
             sex == "female"
           
    )
  ## FILTER sample data for the duplicated measures
  t <- table(sample_data$host_subject_id)
  t <- t[t == 1]
  s <- (sample_data[sample_data$host_subject_id %!in% names(t),])
  s <- s[order(s$host_subject_id),]
  sample_data <- sample_data[-which(rownames(sample_data) %in% rownames(s)),]
  otu_table <- (otu_table(physeq) %>% as("matrix"))[,rownames(sample_data)]
  sample_data <- sample_data[rownames(sample_data) %in% colnames(otu_table),]
  indg1 <- which(sample_data$diagnosis == "no")
  indg2 <- which(sample_data$diagnosis == "CD")
  sample_data$Group <- sample_data$diagnosis
  sample_data <- sample_data[rownames(sample_data) %in% colnames(otu_table),]
  otu_table <- otu_table[rowSums(otu_table)>0,]
  design <- model.matrix(as.formula("~ Group"), data = as.data.frame(as(sample_data,"matrix")))
}
## Test if Demographic Variables May Be Controlled For

# fit <- glm(sample_data$diagnosis ~ sample_data$sex
# family = "binomial") # Alt Model
# coef(summary(fit))

# fit <- glm(sample_data$diagnosis ~ 1, 
# family = "binomial") # Null Model
# coef(summary(fit))

################################################################################

## Normalizations

## Apply rank normalization, analyze distributions. 
ranks <- apply(otu_table, 2, rank)

## TSS
tss <- apply(otu_table, 2, function(x)x/sum(x))

## CSS
library(metagenomeSeq)
aux <- newMRexperiment(counts = as(otu_table,"matrix"))
normFacts <- metagenomeSeq::calcNormFactors(
  obj = aux, p = cumNormStatFast(aux))
w <- drop(as.matrix(normFacts))
css <-  t(t(otu_table) / w)

## TMM
library(edgeR)
normFacts <- edgeR::calcNormFactors(otu_table, method = "TMM")
w = normFacts * colSums(otu_table)
tmm <-  t(t(otu_table) / w)


################################################################################

## Run Differential Abundance Tests

## store rejected OTU by each methodology
rej_list <- list()

## Simple Tests
otu_norms <- list("ranks" = ranks,
                  "tss" = tss,
                  "css" = css,
                  "tmm" = tmm,
                  "none" = otu_table)
for(i in 1:length(otu_norms)){
  otu_norm <- otu_norms[[i]]
  pvals <- fastT(otu_norm, indg1, indg2, sign = F) 
  pvals <- pvals[pvals > 0]
  padj <- p.adjust(abs(pvals), "BH")
  rej <- which(padj < .05)
  rej_list[[paste0("t-test_",names(otu_norms)[i])]] <- padj[names(rej)]
  
  if(names(otu_norms)[i] != "ranks"){
    pvals <- apply(t(otu_norm), 2, function(x){
      wilcox.test(x[indg1], x[indg2],exact = F)$p.value
    })
    
    
    padj <- p.adjust(abs(pvals), "BH")
    rej <- which(padj < .05)
    rej_list[[paste0("wilcox_",names(otu_norms)[i])]] <- padj[names(rej)]
  }
}

## ALDEx2
library(ALDEx2)
conditions <- rep(NA, ncol(otu_table))
conditions[indg1] <- "indg1"
conditions[indg2] <- "indg2"
aldex <- aldex(as.data.frame(round(otu_table)), 
               conditions = conditions,
               test ="t", 
               effect = FALSE)
pvals <- aldex$we.ep
padj <- p.adjust(pvals, "BH")
names(padj) <- rownames(aldex)
rej <- which(padj < 0.05)
rej_list[[paste0("aldex")]] <- padj[names(rej)]

## Limma Voom with TMM
library(limma)
otu_norm <- tmm
res <- voom(counts = otu_norm, design = design, plot = FALSE) %>%
  lmFit(design = design) %>%
  eBayes(robust = T) %>%
  topTable(coef = 2, n = nrow(otu_norm), sort.by = "none")
pval <- res$P.Value
names(pval) <- rownames(res)
padj <- p.adjust(pval, "BH")
rej <- which(padj < 0.05)
rej_list[[paste0("voom_tmm")]] <- padj[names(rej)]

## Limma Voom with TSS
library(limma)
otu_norm <- tss
res <- voom(counts = otu_norm, design = design, plot = FALSE) %>%
  lmFit(design = design) %>%
  eBayes(robust = T) %>%
  topTable(coef = 2, n = nrow(otu_norm), sort.by = "none")
pval <- res$P.Value
names(pval) <- rownames(res)
padj <- p.adjust(pval, "BH")
rej <- which(padj < 0.05)
rej_list[[paste0("voom_tss")]] <- padj[names(rej)]

## Limma Voom with CSS
library(limma)
otu_norm <- css
res <- voom(counts = otu_norm, design = design, plot = FALSE) %>%
  lmFit(design = design) %>%
  eBayes(robust = T) %>%
  topTable(coef = 2, n = nrow(otu_norm), sort.by = "none")
pval <- res$P.Value
names(pval) <- rownames(res)
padj <- p.adjust(pval, "BH")
rej <- which(padj < 0.05)
rej_list[[paste0("voom_css")]] <- padj[names(rej)]

## MetagenomeSeq
library(metagenomeSeq)
otu_norm <- css
MGS <- newMRexperiment(counts = otu_norm, phenoData = AnnotatedDataFrame(sample_data), 
                       normFactors = colSums(otu_norm))
res <- fitZig(MGS, design) %>%
  MRfulltable(number = nrow(get("counts", assayData(MGS))))
pvals <- res[, c("pvalues")]
names(pvals) <- rownames(res)
pvals <- pvals[!is.nan(pvals)]
padj <- p.adjust(pvals, "BH")
rej <- which(padj < 0.05)
rej_list[[paste0("metagenomeSeq")]] <- padj[names(rej)]

## DESeq2 
library(DESeq2)
geoMeans <- apply(round(otu_table), 
                  1, 
                  function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
res <- DESeqDataSetFromMatrix(
  countData = round(otu_table),
  colData = as.data.frame(sample_data$Group),
  design = design) %>%
  estimateSizeFactors(geoMeans = geoMeans) %>%
  DESeq %>%
  results
pvals <- res$pvalue
names(pvals) <- rownames(res)
pvals <- pvals[!is.na(pvals)]
padj <- p.adjust(pvals, "BH")
rej <- which(padj < 0.05)
rej_list[[paste0("DESeq2")]] <- padj[names(rej)]


## edgeR
library(edgeR)
res <- DGEList(counts = otu_table, group = sample_data$Group) %>%
  calcNormFactors(method = "TMM") %>%
  estimateGLMRobustDisp %>%
  glmQLFit(robust = TRUE,design = design) %>%
  glmQLFTest(coef = 2)
pvals <- res$table$PValue
padj <- p.adjust(pvals, "BH")
names(padj) <- rownames(res$table)
rej <- which(padj < 0.05)
rej_list[[paste0("edgeR")]] <- padj[names(rej)]

################################################################################

## Shared Rejections (Significant Taxa) Between Methods

## the proportion of shared OTU between methodologies m1 and m2, 
## out of the total number rejections from methodology m1
intersect_prop <- function(m1,m2){
  n1 <- names(rej_list[[m1]])
  n2 <- names(rej_list[[m2]])
  (intersect(n1,n2) %>% length)/length(n1) %>% round(4)
}

## Intersect Matrix (proportion of each row out of the shared taxa between columns)
shared_taxa <- matrix(c(1:length(rej_list)),nrow = length(rej_list),ncol = length(rej_list))
for(i in 1:length(rej_list)){
  for(j in 1:length(rej_list)){
    shared_taxa[i,j] <- intersect_prop(i,
                                       j)
  }
}
rownames(shared_taxa) <- names(rej_list)
colnames(shared_taxa) <- names(rej_list)
shared_taxa[is.nan(shared_taxa)] <- 0
shared_taxa

################################################################################

## Display Results
otu_norm <- ranks
rej <- rej_list$`t-test_ranks`

## For Zeller
if(dataset == "Zeller"){
  rej_taxa <- names(rej)
  TaxaIdentified <- data.frame()
  for(taxa in names(rej)){
    x <- otu_norm[taxa,indg1]
    y <- otu_norm[taxa,indg2]
    t <- t.test(x,y)
    TaxaIdentified <- TaxaIdentified %>% rbind(c(
      t$estimate[c(2,1)],
      t$estimate[1] - t$estimate[2],
      t$conf.int*(sign(t$estimate[1] - t$estimate[2])),
      t$p.value,
      rej[taxa])) %>%
      as.data.frame %>%
      set_colnames(c("Mean Rank Control",
                     "Mean Rank CRC",
                     "Mean Difference",
                     "LB (2.5%)",
                     "UB (97.5%)",
                     "Pvalue",
                     "Padj"))
  }
  TaxaIdentified$taxon<- names(rej)
  TaxaIdentified <- TaxaIdentified[,c(ncol(TaxaIdentified),1:(ncol(TaxaIdentified)-1))]
  TaxaIdentified <- TaxaIdentified[order(TaxaIdentified$Padj),]
  options(scipen = 4)
  View(TaxaIdentified)
}


## For RISK
if(dataset == "RISK"){
  rej <- rej_list$`t-test_ranks`
  rej_taxa <- tax[names(rej),] %>% as.data.frame %>% mutate(padj = rej)
  TaxaIdentified <- data.frame()
  for(taxa in names(rej)){
    x <- otu_norm[taxa,indg1]
    y <- otu_norm[taxa,indg2]
    t <- t.test(x,y)
    TaxaIdentified <- TaxaIdentified %>% rbind(c(
      t$estimate[c(1,2)],
      t$estimate[2] - t$estimate[1],
      -rev(t$conf.int),
      # t$conf.int*(sign(t$estimate[2] - t$estimate[1])),
      t$p.value,
      rej[taxa])) %>%
      as.data.frame %>%
      set_colnames(c("Mean Rank Control",
                     "Mean Rank CD",
                     "Mean Difference",
                     "LB (2.5%)",
                     "UB (9.75%)",
                     "Pvalue",
                     "Padj"))
  }
  TaxaIdentified$taxon<- tax[names(rej),c("Order","Family","Genus","Species")] %>%
    apply(1,function(x){
      x[is.na(x)] <- "Unknown"
      paste(x,collapse = " ")
    })
  TaxaIdentified <- TaxaIdentified[,c(ncol(TaxaIdentified),1:(ncol(TaxaIdentified)-1))]
  TaxaIdentified <- TaxaIdentified[order(TaxaIdentified$Padj),]
  options(scipen = 4)
  View(TaxaIdentified)
}

################################################################################ 

## Compare to Original Findings

## Compare Zeller to Original Findings
if (dataset == "Zeller") {
  original_findings <- c(            # explicitly mentioned taxa from Figure 1A
    "unclassified Fusobacterium",    # found by ranks
    "unclassified Fusobacterium",    # found by ranks
    "Fusobacterium nucleatum",       # found by ranks
    "Peptostreptococcus stomatis",
    "Porphyromonas asaccharolytica", # found by ranks
    "Clostridium symbiosum",
    "Clostridium hylemonae",
    "Bacteroides fragilis",
    "Lactobacillus salivarius",
    "Fusobacterium gonidiaformans",
    "Lactobacillus ruminis",
    "Eubacterium rectale",
    "Bacteroides caccae",
    "Eubacterium ventriosum",
    "Clostridium scindens",
    "Eubacterium eligens",
    "Bifidobacterium angulatum",
    "Methanosphaera stadtmanae",
    "Dorea formicigenerans",
    "Butyrivibrio crossotus",
    "Phascolarctobacterium succinatutens",
    "unnamed Ruminococcus",
    "Streptococcus salivarius"
  )
  
  ## Match to Taxa
  rownames <- sapply(rownames(otu_table), function(x) {
    paste(strsplit(x, " ")[[1]][1:2], collapse = " ")
  })
  species <- rownames[which(rownames %in% original_findings)]
  
  ## Rough Estimates of Performance
  RoughRes <- data.frame(
    "Percent Agreement" = lapply(rej_list, function(x)
      mean(names(x) %in% names(species))) %>% unlist %>% sapply(function(x)
        ifelse(is.nan(x), 0, x)) %>% unlist,
    "Number Rejected" = lapply(rej_list, length) %>% unlist,
    "Percent Shared with T-Test Ranks" = lapply(rej_list, function(x) 
      mean(names(x) %in% (names(rej_list$`t-test_ranks`)))) %>% unlist %>% sapply(function(x)
        ifelse(is.nan(x), 0, x)) %>% unlist
  )
  RoughRes <- RoughRes[order(rownames(RoughRes),decreasing = T), ]
  print(RoughRes[order(rownames(RoughRes),decreasing = T), ])
}

## For RISK
if (dataset == "RISK") {
  
  # explicitly mentioned taxa from supplementary File 1 PDF
  increased <- c(
    "Enterobacteriaceae", # found by ranks
    "Pasteurellaceae",    # found by ranks
    "Fusobacteriaceae",   
    "Neisseriaceae",      
    "Veillonellaceae",    # found by ranks
    "Gemellaceae")        
  
  decreased <- c(
    "Bacteroidales",       # found by ranks
    "Clostridiales",       # (excluding Veillonellaceae) found by ranks
    "Erysipelotrichaceae", # found by ranks
    "Bifidobacteriaceae")
  
  ## Recognized Taxa
  original_findings <-
    rownames(tax[apply(tax, 1, function(x)
      any(as.character(x) %in% c(increased, decreased))), ])
  
  ## Rough Estimates of Performance
  RoughRes <- data.frame(
    ## This is the estimate of 1-FDR
    "Percent Agreement" = lapply(rej_list, function(x) 
      mean(names(x) %in% (original_findings))) %>% unlist %>% sapply(function(x)
        ifelse(is.nan(x), 0, x)) %>% unlist,
    "Number Rejected" = lapply(rej_list, length) %>% unlist,
    "Percent Shared with T-Test Ranks" = lapply(rej_list, function(x) 
      mean(names(x) %in% (names(rej_list$`t-test_ranks`)))) %>% unlist %>% sapply(function(x)
        ifelse(is.nan(x), 0, x)) %>% unlist
  )
  RoughRes <- RoughRes[order(rownames(RoughRes),decreasing = T), ]
  print(RoughRes)
  
  ## All of these taxa were listed by original findings except Turibacter
  verify_ind <- sapply(c(increased, decreased), function(x) {
    sapply(x, function(y) {
      c(
        which(rej_taxa$Family == y),
        which(rej_taxa$Genus == y),
        which(rej_taxa$original_findings == y),
        which(rej_taxa$Order == y)
      )
    }) %>% unlist %>% unique
  }) %>% unlist %>% unique
  rej_taxa_rem <- rej_taxa[-verify_ind, ]
  View(rej_taxa_rem)
  
  ## Greater Turicibacter is a novel finding
  t <- table(rej_taxa_rem$Order)
  t <- t[t > 0]
  TaxaIdentified[sapply(names(t[order(t)]), function(x)
    grep(x, TaxaIdentified$taxon)) %>% unlist, ]
}


## Record/Write Results
options(digits = 2)
options(scipen = 2)
RoughRes$Percent.Agreement <- RoughRes$Percent.Agreement * 100
RoughRes$Percent.Shared.with.T.Test.Ranks <- RoughRes$Percent.Shared.with.T.Test.Ranks * 100
write.csv(RoughRes, file = paste0(dataset, "_RoughRes.csv"))

TaxaIdentified <- TaxaIdentified[order(TaxaIdentified$`Mean Difference`,decreasing = TRUE),]
View(TaxaIdentified)
write.csv(TaxaIdentified,file = paste0(dataset,"_Rejects.csv"))
################################################################################

## Nonparametric Permutation T-Test ##
permTest <- function(    otu_table,    # rows as taxa, columns as samples
                         ncores = 8, # number of cores used
                         N = 1000,   # number of permutations
                         indg1,      # column indices of group 1
                         indg2,      # column indices of group 2
                         alpha = 0.05, # significance level 
                         adj = "",   # BH? FDR? what kind of adjustment, if any
                         makeCluster = T, # if we should make a cluster within the function
                         juststat = FALSE){ # return only the permuted t-statistics
  ## Observed T-Statistic
  l1 <- length(indg1)
  l2 <- length(indg2)
  S1 <- (rowSums((otu_table[,indg1] - rowMeans(otu_table[,indg1]))^2) / (l1 - 1))
  S2 <- (rowSums((otu_table[,indg2] - rowMeans(otu_table[,indg2]))^2) / (l2 - 1))
  sigmas <- sqrt ( S1/l1 + S2/l2 )
  TObs <- (rowMeans(otu_table[,indg1]) - rowMeans(otu_table[,indg2])) / sigmas
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
    S1 <- (rowSums((otu_tablePerm[,indg1] - rowMeans(otu_tablePerm[,indg1]))^2) / (l1 - 1))
    S2 <- (rowSums((otu_tablePerm[,indg2] - rowMeans(otu_tablePerm[,indg2]))^2) / (l2 - 1))
    sigmas <- sqrt ( S1/l1 + S2/l2 )
    TP <- (rowMeans(otu_tablePerm[,indg1]) - rowMeans(otu_tablePerm[,indg2])) / sigmas
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

## Run Permutation t-test non-parametric alternative to Welchs t-test
library(foreach)
library(doSNOW)
library(parallel)

cl2 <- makeCluster(4)
registerDoSNOW(cl2)
rej <- permTest(ranks, 
                N = 500,
                indg1 = indg1,
                indg2 = indg2,
                adj = "BH")
stopCluster(cl2)
rej


