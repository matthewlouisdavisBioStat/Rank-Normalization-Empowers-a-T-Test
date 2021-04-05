# Rank Normalization Empowers a t-test for Microbiome Differential Abundance Analysis while Controlling for False Discoveries.
## Supplementary Materials and Software Associated with the Publication 

Attached here are all of the data files and code needed to reproduce all the findings and results from the study. Contact matthew-l-davis@uiowa.edu for any other requests.

https://doi.org/10.1093/bib/bbab059

##  SimpleExample.R
  
Performs rank normalization with a t-test on fabricated microbiome data, only using base R.

##  RCode

Much of the following code was adapted from materials previously made freely available by Hawinkel et al. at https://users.ugent.be/~shawinke/ABrokenPromise/. See their original benchmarking study for further details https://pubmed.ncbi.nlm.nih.gov/28968702/.

- DataGeneration.R: Generate simulated datasets under 'H1 without compensation' using Mid Vagina HMP template, with negative binomial and beta-binomial distributions.
- ImportHMPData.R: Import/Clean HMP data for simulations.
- Normalizations.R: Code to perform normalizations for simulation testing,called by RunTestFunction.R.
- PermutationTestExample.R: Code to perform permutation t-test in parallel on fabricated dataset (as a non-parametric alternative to Welch's t-test).
- PlotResults1point5.R: Plot and recreate the scatterplots summarizing simulation results (with fold-change of 1.5).
- PlotResultsLargeSmall.R: Plot and recreate the scatterplots summarizing simulation results (with fold-changes of 3 and 5).
- RealDataAnalysis.R: Code to perform all real data analysis for all methodologies.
- RunSim.R: Run differential abundance tests on simulated datasets, and record performance.
- RunTestFunction.R: Called by RunSim.R, to apply the differential abundance tests to a given dataset.
- Tests.R: Differential abundance detection methodologies for simulations, called by RunTestFunction.R.
- fastTTest: Perform a t-test simultaneously across all rows of an OTU table.
- msWALDHMP.R: Code copied and pasted from https://github.com/mafed/msWaldHMP for modelling microbiome counts, called by DataGeneration.R.
- RankNormPlots.Rmd: RMarkdown code to plot and explore the ranks and t-statistics on real datasets
- rankTests (for running sim).R: Called by RunTestFunction.R to perform rank normalization and a t-test on simulated datasets.
- permTests (for running sim).R: (Optionally) called by RunTestFunction.R to perform rank normalization and a permutation t-test on simulated datasets.

##  data

Data Sources: 
RISK (as curated by MicrobeDS https://github.com/twbattaglia/MicrobeDS)
Zeller (see https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4299606/) 
HMP (see https://rdrr.io/github/waldronlab/HMP16SData/src/inst/scripts/make-data.R)

- covList.RData: Variance-covariance structure of Mid Vigina HMP Template as estimated by SpiecEasi (loaded in DataGeneration.R). This was requested directly from Hawinkel, who kindly provided it.
- msb0010-0766-sd2.xlsx: Sample Data/Metadata on participants in Zeller study (loaded by RealDataAnalysis.R.).
- msb0010-0766-sd4.xlsx: OTU table in Zeller study, taken from supplementary materials (loaded by RealDataAnalysis.R).
- otu_table_psn_v13.txt.gz and otu_table_psn_v35.txt.gz: OTU tables for HMP data (loaded by ImportHMPData.R).
- otu_table_psn_v13.txt.gz and otu_table_psn_v35.txt.gz: OTU tables for HMP data (loaded by ImportHMPData.R.).
- v13_map_uniquebyPSN.txt.bz and v35_map_uniquebyPSN.txt.bz: Connects metadata/sample data for HMP data to OTU table (loaded by ImportHMPData.R). 
- v13_v35_psn_otu.genus.fixed.zip: zipped files of OTU tables for HMP data (both v13_psn_otu.genus.fixed.txt and v13_psn_otu.genus.fixed.txt loaded by ImportHMPData.R).

##  results

Data tables corresponding to simulation and real data analysis results presented in the manuscript. 

- MainSimulationResults.csv: tables of simulation performance for all differential abundance detection methodologies for all simulations conducted. Data here mirrors what is presented by the simulation results scatterplots. 
- MainSimulationResultsLarge_NegBin.pdf and MainSimulationResultsLarge_BetaBin.pdf: Figures presented in the simulation results, showing the Sensitivity vs. 1-FDR scatter plots of differential abundance detection methodologies, performed on simulated datasets with the marginal distritions of OTU following a zero-inflated negative binomial distribution and beta-binomial distribution respectively at sample sizes of 50, 100, and 150 with fold-changes of 3 and 5.
- MainSimulationResultsSmall_NegBin.pdf and MainSimulationResultsSmall_BetaBin.pdf: Figures presented in the simulation results, showing the Sensitivity vs. 1-FDR scatter plots of differential abundance detection methodologies, performed on simulated datasets with the marginal distritions of OTU following a zero-inflated negative binomial distribution and beta-binomial distribution respectively at sample sizes of 5, 15, and 25 with fold-changes of 3 and 5.
- S1_MainSimulationResults1.5_BetaBin andS S2_MainSimulationResults1.5_NegBin: Figures presented as supplementary materials, showing the Sensitivity vs. 1-FDR scatter plots of differential abundance detection methodologies, performed on simulated datasets with the marginal distritions of OTU following a zero-inflated negative binomial distribution and beta-binomial distribution respectively at sample sizes of 5, 15, 25, 50, 100, and 150 with fold-change of 1.5. These are referred to as Supplementary Figures "S1" and "S2" in the manuscript. 
- RISK_Rejects.csv and Zeller_Rejects.csv: csv files containing the OTU identified as significantly differentially abundant by rank normalization with a t-test on RISK and Zeller data respectively, along with estimated mean differences in ranks, 95% confidence intervals for this mean difference, and raw as well as BH-adjusted p-values.
- RISK_RoughRes.csv and Zeller_RoughRes.csv: csv files containing results from real data analysis on RISK and Zeller data respectively, comparing for each differential abundance detection methodology the % of OTU rejected that were identified by the original findings, the number of OTU rejected (N), and the % of OTU rejected that were also rejected by a t-test paired with rank normalization. 
- Ranks.png: Histograms of observed ranks for 3 microbiome datasets as presented in the Methods section of the manuscript.
