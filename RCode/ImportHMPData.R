
  ## I'm using R 3.5.1 Right Now

## All code is adapted from materials made available from
##     https://users.ugent.be/~shawinke/ABrokenPromise/

## Get data files from
# OTU table:
# http://downloads.ihmpdcc.org/data/HMQCP/otu_table_psn_v13.txt.gz
# http://downloads.ihmpdcc.org/data/HMQCP/otu_table_psn_v35.txt.gz

# Sample metadata:
# http://downloads.ihmpdcc.org/data/HMQCP/v13_map_uniquebyPSN.txt.bz2
# http://downloads.ihmpdcc.org/data/HMQCP/v35_map_uniquebyPSN.txt.bz2

# Subject metadata files with SRS identifiers:
# https://www.hmpdacc.org/hmp/doc/ppAll_V13_map.txt
# https://www.hmpdacc.org/hmp/doc/ppAll_V35_map.txt

# Phylogeny:
# http://downloads.ihmpdcc.org/data/HMQCP/rep_set_v13.tre.gz
# http://downloads.ihmpdcc.org/data/HMQCP/rep_set_v35.tre.gz



###################################################################################################

  ## Working directory, load packages
WD =  "C:/Users/Matthew/Documents/Courses/Kai/Final Results/Final Plots/DataGeneration"
knitr::opts_chunk$set(cache = TRUE, autodep = TRUE, root.dir = WD,tidy=TRUE, cache.lazy=FALSE)
setwd(WD)
dataWD <- "C:/Users/Matthew/Documents/Courses/Kai/Final Results/Final Plots/DataGeneration/"
# The required package list:
reqpkg = c("phyloseq","reshape2")
# Load all required packages and show version
for(i in reqpkg)
{
  print(i) 
  print(packageVersion(i))
  library(i, quietly=TRUE, verbose=FALSE, warn.conflicts=FALSE, character.only=TRUE)
}

## Load metadata
sampleDataV13 <- read.delim(paste0(dataWD, "v13_map_uniquebyPSN.txt.bz2"))
sampleDataV35 <- read.delim(paste0(dataWD, "v35_map_uniquebyPSN.txt.bz2"))

# # setbody region
regionFun = function(df) {
  within(df, {
    id_vagina = HMPbodysubsite %in% c("Mid_vagina", "Posterior_fornix", 
                                      "Vaginal_introitus")
    id_oralCavity = HMPbodysubsite %in% c("Throat", "Hard_palate", "Keratinized_gingiva", 
                                          "Saliva", "Subgingival_plaque", "Supragingival_plaque", "Palatine_Tonsils", 
                                          "Tongue_dorsum", "Buccal_mucosa")
    id_skin = HMPbodysubsite %in% c("L_Retroauricular_crease", "R_Retroauricular_crease", 
                                    "L_Antecubital_fossa", "R_Antecubital_fossa")
    id_airways = HMPbodysubsite %in% c("Anterior_nares")
    id_blood = HMPbodysubsite %in% c("Blood")
    id_stool = HMPbodysubsite %in% c("Stool")
    HMPregion = factor("character", levels = c("Vagina", "Oral Cavity", 
                                               "Skin", "Airways", "Blood", "Gut"))
    HMPregion[id_vagina] = "Vagina"
    HMPregion[id_oralCavity] = "Oral Cavity"
    HMPregion[id_skin] = "Skin"
    HMPregion[id_airways] = "Airways"
    HMPregion[id_blood] = "Blood"
    HMPregion[id_stool] = "Gut"
    HMPregion = as.factor(HMPregion)
    id_vagina = id_oralCavity = id_skin = id_airways = id_blood = id_stool = NULL
  })
}

## Prepare metadata
sampleDataV13 = regionFun(sampleDataV13)
sampleDataV35 = regionFun(sampleDataV35)
rownames(sampleDataV13) <- sampleDataV13$SampleID
rownames(sampleDataV35) <- sampleDataV35$SampleID
metaVars2Keep <- c("HMPbodysubsite", "HMPregion", "sex", "RSID", "RUNCENTER")
sampleDataV13 <- sampleDataV13[, metaVars2Keep]
sampleDataV35 <- sampleDataV35[, metaVars2Keep]
sampleDataV13 <- sample_data(sampleDataV13)
sampleDataV35 <- sample_data(sampleDataV35)
rankNames <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
sitesV13 <- names(table(sampleDataV13$HMPbodysubsite))
sitesV35 <- names(table(sampleDataV35$HMPbodysubsite))
if (!all.equal(sitesV13, sitesV35)) {
  stop("Mismatch in sample data coming from the two regions V1-3 and V3-5!")
} else {
  invisible()
}


## Load OTU data
if (!file.exists(paste0(dataWD,"physeqListV13.RData"))) {
  # # # OTU table for region V13
  otuTabV13 <- read.table(paste0(dataWD, "v13_psn_otu.genus.fixed.txt"), header = TRUE, 
                          sep = "\t", row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE, 
                          comment.char = "#")
  # colClasses = c('character', rep.int('integer', 2910), 'character'))
  
  colnames(otuTabV13)[ncol(otuTabV13)] <- "Consensus.Linage"
  linageV13 <- gsub(pattern = ";", replacement = "|", x = otuTabV13$Consensus.Linage, 
                    fixed = TRUE)
  otuTabV13 <- otuTabV13[, colnames(otuTabV13) != "Consensus.Linage"]
  
  indOTUs2Keep <- rowSums(otuTabV13) > 0L
  otuTabV13 <- otuTabV13[indOTUs2Keep, ]
  linageV13 <- linageV13[indOTUs2Keep]
} else {
  invisible()
}

## Tax Table
if (!file.exists(paste0(dataWD,"physeqListV13.RData"))) {
  taxMatrixV13 <- as.matrix(colsplit(string = linageV13, pattern = "\\|", 
                                     names = rankNames))
  rownames(taxMatrixV13) <- rownames(otuTabV13)
} else {
  invisible()
}


## Ensure rownames of metadata and colnames of OTU data match
rownames(sampleDataV13) <- substr(rownames(sampleDataV13), 3, nchar(rownames(sampleDataV13)))
colnames(otuTabV13) <- substr(colnames(otuTabV13), 3, nchar(colnames(otuTabV13)))

##
if (!file.exists(paste0(dataWD,"physeqListV13.RData"))){
  physeqListV13 <- lapply(sitesV13, FUN = function(siteRun) {
    sampleInd <- rownames(sampleDataV13)[sampleDataV13$HMPbodysubsite == 
                                           siteRun]
    sampleInd <- colnames(otuTabV13)[colnames(otuTabV13) %in% sampleInd]
    # tmpOtu <- simpleTrim(otuTabV13[, sampleInd])
    tmpOtu <- otuTabV13[, sampleInd]
    
    ## trimming by 1% prevalence on 1 counts
    prevalence <- rowMeans(tmpOtu > 1L)
    tmpOtu <- tmpOtu[prevalence > 0.01, ]
    ## trim libraries with size < 100
    tmpOtu = tmpOtu[, colSums(tmpOtu) >= 100]
    phyloseq(otu_table(tmpOtu, taxa_are_rows = TRUE), 
             sample_data(sampleDataV13[sampleInd,]),
             tax_table(taxMatrixV13[rownames(otuTabV13),]))
  })
  names(physeqListV13) <- gsub(pattern = "_", replacement = ".", x = sitesV13, 
                               fixed = TRUE)
  # utils:::print.object_size(object.size(physeqListV13), unit = 'MB')
  rm(otuTabV13, sampleDataV13, taxMatrixV13)
  gc(reset = TRUE)
  save(physeqListV13, file = paste0(dataWD, "physeqListV13.RData"))
} else {
  invisible()
}


## Repeat the above steps in one block of code for V35 region
if (!file.exists(paste0(dataWD, "physeqListV35.RData"))) {
  # # # OTU table for region V35
  otuTabV35 <- read.table(paste0(dataWD, "v35_psn_otu.genus.fixed.txt"), header = TRUE, 
                          sep = "\t", row.names = 1L, check.names = FALSE, stringsAsFactors = FALSE, 
                          comment.char = "#")
  # colClasses = c('character', rep.int('integer', 2910), 'character'))
  
  colnames(otuTabV35)[ncol(otuTabV35)] <- "Consensus.Linage"
  linageV35 <- gsub(pattern = ";", replacement = "|", x = otuTabV35$Consensus.Linage, 
                    fixed = TRUE)
  otuTabV35 <- otuTabV35[, colnames(otuTabV35) != "Consensus.Linage"]
  
  indOTUs2Keep <- rowSums(otuTabV35) > 0L
  otuTabV35 <- otuTabV35[indOTUs2Keep, ]
  linageV35 <- linageV35[indOTUs2Keep]
  #utils:::print.object_size(object.size(otuTabV35), unit = "MB")
  taxMatrixV35 <- as.matrix(colsplit(string = linageV35, pattern = "\\|", 
                                     names = rankNames))
  rownames(taxMatrixV35) <- rownames(otuTabV35)
  rownames(sampleDataV35) <- substr(rownames(sampleDataV35), 3, nchar(rownames(sampleDataV35)))
  colnames(otuTabV35) <- substr(colnames(otuTabV35), 3, nchar(colnames(otuTabV35)))
  physeqListV35 <- lapply(sitesV35, FUN = function(siteRun) {
    sampleInd <- rownames(sampleDataV35)[sampleDataV35$HMPbodysubsite == 
                                           siteRun]
    sampleInd <- colnames(otuTabV35)[colnames(otuTabV35) %in% sampleInd]
    # tmpOtu <- simpleTrim(otuTabV35[, sampleInd])
    tmpOtu <- otuTabV35[, sampleInd]
    ## trimming by 1% prevalence on 1 counts
    prevalence <- rowMeans(tmpOtu > 1L)
    tmpOtu <- tmpOtu[prevalence > 0.01, ]
    ## trim libraries with size < 100
    tmpOtu = tmpOtu[, colSums(tmpOtu) >= 100]
    phyloseq(otu_table(tmpOtu, taxa_are_rows = TRUE), sample_data(sampleDataV35[sampleInd, 
                                                                                ]), tax_table(taxMatrixV35[rownames(tmpOtu), ]))
  })
  names(physeqListV35) <- gsub(pattern = "_", replacement = ".", x = sitesV35, 
                               fixed = TRUE)
  utils:::print.object_size(object.size(physeqListV35), unit = "MB")
  rm(otuTabV35, sampleDataV35, taxMatrixV35)
  gc(reset = TRUE)
  save(physeqListV35, file = paste0(dataWD,"physeqListV35.RData"))
} else {
  invisible()
}





