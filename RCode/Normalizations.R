
    ## Functions for Normalization
    ## code is adapted from from https://users.ugent.be/~shawinke/ABrokenPromise/03_diffAbundDetect.html


### Rarefy
normRare <- function(physeq)
{
  physeq <- rarefy_even_depth(physeq, trimOTUs = FALSE)
  aux <- data.frame(sample_data(physeq))
  aux$"NF.rare" <- 1
  sample_data <- sample_data(aux)
  rownames(sample_data) <- rownames(sample_data(physeq))
  sample_data(physeq) <- sample_data
  physeq
}# END 

#- function: normNone
normNone <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.none" <- 1
  sample_data <- sample_data(aux)
  rownames(sample_data) <- rownames(sample_data(physeq))
  sample_data(physeq) <- sample_data
  physeq
}# END - function: normNone

### Total sum scaling, also known as proportion normalization, dividing by library sizes
normTSS <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.TSS" <- sample_sums(physeq)#/exp(mean(log(sample_sums(physeq))))
  sample_data <- sample_data(aux)
  rownames(sample_data) <- rownames(sample_data(physeq))
  sample_data(physeq) <- sample_data
  physeq
  
}# END - function: normNone

###Rarefy
normRare <- function(physeq)
{
  physeq <- rarefy_even_depth(physeq)
  aux <- data.frame(sample_data(physeq))
  aux$"NF.rare" <- 1
  sample_data <- sample_data(aux)
  rownames(sample_data) <- rownames(sample_data(physeq))
  sample_data(physeq) <- sample_data
  physeq
}# END - function: normNone

### edgeR normalisations: TMM and RLE
normEdgeR <- function(physeq, method = c('TMM', 'RLE', 'upperquartile'))
{
  # require(edgeR)
  otuTab <- as(otu_table(physeq), "matrix")
  if (!taxa_are_rows(physeq))
  {
    otuTab <- t(otuTab)
  } else {}
  
  if (method == "upperquartile")
  {
    scaledCounts <- t(otuTab) / colSums(otuTab)
    tmpNF <- apply(scaledCounts, MARGIN = 1L, FUN = function(x)
      quantile(x[x != 0], probs = .75))
    normFacts <- tmpNF/exp(mean(log(tmpNF)))
    method <- "UQ"
  } else
  {
    normFacts <- edgeR:::calcNormFactors(otuTab, method = method)
    
  }# END - ifelse: upperquartile only of non-zero counts
  #VERY IMPORTANT: multiply by library sizes and renormalize. edgeR calculates scaling factors, which still have to be multiplied by library sizes to get to the size factors of effective sequencing depth, i.e. robust estimates of the library sizes
  #normFacts = normFacts*sample_sums(physeq)
  #normFacts=normFacts/exp(mean(log(normFacts)))
  if (all(is.na(normFacts))) #Resort to proportion normalization in case of failure for all samples
  {
    normFacts = sample_sums(physeq)
  }
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- paste("NF", method, sep = ".")
  physeq@sam_data@names <- aux
  physeq
}# END - function: normEdgeR


### function that apply different normalisations and build *DESeqDataSet* object
### function that apply different normalisations and build *DESeqDataSet* object
### for DESeq2 analysis
normDESeq2 <- function(physeq, whichOTUs = NULL, method = "ratio")
{
  # require(DESeq2)
  method <- match.arg(method)
  
  ### Coerce count data to vanilla matrix of integers and check if there are zeroes
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  ## select which OTUs to analyse
  if (!missing(whichOTUs) || !is.null(whichOTUs))
  {
    physeq <- prune_taxa(taxa_names(physeq)[whichOTUs], physeq)
  } else {}# END - if: whichOTUs
  
  otuTab <- as(otu_table(physeq), "matrix")
  if (any(otuTab == 0))
  {
    otuTab <- otuTab + 1L
  } else {}
  
  #   otu_table(physeq) <- otu_table(otuTab, taxa_are_rows = TRUE)
  
  ## Calculate size factors
  if (method == "ratio")
  {
    normFacts <- DESeq2::estimateSizeFactorsForMatrix(otuTab)
  } else
  {
    #     normFacts <- DESeq2::estimateSizeFactorsIterate(otuTab)
  }
  
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- "NF.ratio"
  physeq@sam_data@names <- aux
  physeq
}# END - function: normDESeq2

### Cumulative Sum Scaling from *metagenomeSeq*
normCSS <- function(physeq, geoMean = FALSE, rel = 0.1)
{
  # require(metagenomeSeq)
  ### force orientation OTUs x samples
  if (!taxa_are_rows(physeq))
  {
    physeq <- t(physeq)
  } else {}
  
  otuTab <- as(otu_table(physeq), "matrix")
  
  aux <- newMRexperiment(counts = otuTab)
  normFacts <- metagenomeSeq::calcNormFactors(
    obj = aux, p = cumNormStatFast(aux, rel = rel))
  normFacts <- drop(as.matrix(normFacts))
  if (geoMean)
  {
    normFacts <- normFacts / exp(mean(log(normFacts[normFacts > 0]), na.rm = TRUE))
  } else {}
  
  #   physeq@otu_table@.Data <- countsMat * normFac
  
  physeq@sam_data@.Data <- c(physeq@sam_data@.Data, list(normFacts))
  aux <- physeq@sam_data@names
  aux[length(aux)] <- "NF.CSS"
  physeq@sam_data@names <- aux
  return(physeq)
}
                   
 ## Unique for quantile normalization, leave at 1 for now and manually apply in diff abundance test
normQuantile <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.quantile" <- 1
  sample_data(physeq) <- aux
  physeq
}# END - function: normNone
                
                   
  ## This simply labels!! DOES NOT ACTUALLY RANK
normRank <- function(physeq)
{
  aux <- data.frame(sample_data(physeq))
  aux$"NF.rank" <- 1
  sample_data(physeq) <- aux
  physeq
}
