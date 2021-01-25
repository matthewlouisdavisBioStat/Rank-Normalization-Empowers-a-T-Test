
runTests <- function(physeq, truede = truede, degenes = degenes,use_perm = FALSE) {
  returnList = list()
  returnList = within(returnList, {
    
    ## all normalisations
    # 
     physeq <- normEdgeR(physeq = physeq, method = "TMM")
     physeq <- normDESeq2(physeq = physeq)
     physeq <- normCSS(physeq = physeq)
     physeq <- normTSS(physeq = physeq)
     physeq <- normNone(physeq)
    sample_data <- sample_data(physeq)
    sample_data$group <- sample_data$condition
    sample_data(physeq) <- sample_data(sample_data)
    cat("\nNormalisations\t")
    
    # When we rank, we want to normalize before trimming, not after 
    # instead, we will FIRST rank, then only analyze our same trimmming scheme
    source("rankTests (for running sim).R")
    tTest_rank = tTest_rank
    # if(use_perm){
    # source("permTests (for running sim).R")
    # permTest_rank = permTest_rank
    # }
    
    cat("\nrankTest\t")
    # 
     tTest_TSS = applySimpleTests(physeq, test = "t-test", normFacts = "TSS")
     tTest_TMM = applySimpleTests(physeq, test = "t-test", normFacts = "TMM")
     tTest_CSS = applySimpleTests(physeq, test = "t-test", normFacts = "CSS")
     tTest_None = applySimpleTests(physeq, test = "t-test", normFacts = "none")
     wTest_TSS = applySimpleTests(physeq, test = "wilcox", normFacts = "TSS")
     wTest_TMM = applySimpleTests(physeq, test = "wilcox", normFacts = "TMM")
     wTest_CSS = try(applySimpleTests(physeq, test = "wilcox", normFacts = "CSS"), silent = TRUE)
     wTest_None = applySimpleTests(physeq, test = "wilcox", normFacts = "none")
     cat("\nSimple tests\t")
    # 
    # ## limma-voom robust version
    voomTest_TMM = limmaVoomRobust(physeq, normFacts = "TMM")
    cat("\nLimma-Voom robust tests\t")
    # 
    # ## we can extract both t-test and wilcox test from same run
    aldexTTest_none <- aldexTTest(aldexTest(physeq))
     cat("\nALDEx2\t")
    # 
    # ## edgeR robust version
    edgeR_TMM = edgeRRobust(physeq, normFacts = "TMM")

    # ## DESeq2 Default Pipeline
     DESeq2_ratio = try(DESeq2(physeq, normFacts = "ratio"), silent = TRUE)
     if(class(DESeq2_ratio) == "try-error"){
       DESeq2_ratio = try(DESeq2(physeq, normFacts = "ratio"), silent = TRUE)
     }
    cat("\nNB DESeq2 tests\t")
    # 
    # ## metagenomeSeq Zero-Inflated Gaussian
    mgsZig_CSS <- try(metagenomeSeqZIG(physeq, normFacts = "CSS"), silent = TRUE)
    if(class(mgsZig_CSS) == "try-error"){
      mgsZig_CSS <- try(metagenomeSeqZIG(physeq, normFacts = "CSS"), silent = TRUE)
    }
    rm(physeq)
    rm(sample_data)
  })
  cat("Adj. WMW-PI tests\t")
  return(returnList)
}