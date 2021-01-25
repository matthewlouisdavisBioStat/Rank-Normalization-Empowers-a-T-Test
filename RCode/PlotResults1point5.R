
'%!in%' <- function(x, y)
  ! (x %in% y)
## Run Tests
setwd(
  "C:/Users/Matthew/Documents/Courses/Kai/Final Results/Final Plots/DifferentialAbundanceTests"
)
datawd <- paste0(getwd(),"/")
## Load simulation results
load("C:/Users/Matthew/Documents/Courses/Kai/Final Results/Final Plots/DifferentialAbundanceTests/MasterResBrokenPromiseDataAll1point5sBetaBin.RData")
MasterRes$sampleType <- "HMP Mid Vagina"
distr_label <- "BetaBin"


library(magrittr)
library(phyloseq)
library(ALDEx2)
library(limma)
library(edgeR)
library(DESeq2)
library(metagenomeSeq)
library(reshape2)
library(ggplot2)
library(ggpubr)
#r1 is sd x direction
#r2 is sd y direction
geom_bubble <-
  function(r1,
           r2,
           xc,
           yc,
           color = "black",
           fill = NA,
           ...) {
    x <-
      rowMaxs(cbind(rowMins(cbind(
        xc + r1 * cos(seq(0, pi, length.out = 100)), 1
      )), 0))
    ymax <-
      rowMaxs(cbind(rowMins(cbind(
        yc + r2 * sin(seq(0, pi, length.out = 100)), 1
      )), 0))
    ymin <-
      rowMaxs(cbind(rowMins(cbind(
        yc + r2 * sin(seq(0,-pi, length.out = 100)), 1
      )), 0))
    annotate(
      "ribbon",
      x = x,
      ymin = ymin,
      ymax = ymax,
      color = color,
      fill = fill,
      ...
    )
  }

################################################################################

MasterRes$sim <- "No Compensation"
incl_test_norms <- c(
  "t-test_rank",
  #"perm-test_rank",
  "Voom_tmm",
  "t-test_tss",
  #"t-test_none",
  "t-test_css",
  "t-test_tmm",
  "wilcox_css",
  #"wilcox_none",
  "wilcox_tss",
  "wilcox_tmm",
  "ALDEx2",
  "edgeR_tmm",
  "DESeq2_rle",
  "metagenomeSeq_css"
)


tests <-
  c(
    "wilcox",
    "t-test",
    #"perm-test",
    "Voom",
    "rankUnifTest",
    "edgeR",
    "metagenomeSeq",
    "DESeq2",
    "wilcoxUnifTest"
  )
norms <- c("rank", "none")
test_norms <- c(
  #"perm-test_rank",
  "t-test_rank-simple",
  "edgeR_tmm",
  "DESeq2_ratio",
  "metagenomeSeq_css",
  "aldexT_aldex",
  "aldexW_aldex",
  "Voom_tmm",
  "wilcox_css",
  "wilcox_none"
)
test_norms <- unique(MasterRes$test_norm)
xvar <- "test_norm"
yvar <- "sens"
fillvar <- "test_norm"
expandgrid1 <- "m"
expandgrid2 <- "effect"
facet <- "sim"
ymaxlim <- .4

samplesize_bins <-
  list(
    small = c(min(MasterRes$m), 10),
    med = c(15, 25),
    large = c(26, max(MasterRes$m))
  )
effectsize_bins <-
  list(small = c(min(MasterRes$effect), 3), large = c(4, max(MasterRes$effect)))
propde_bins <-
  list(
    small = c(0, .05),
    med = c(.051, .1),
    large = c(.101, 1)
  )
samplesize_bins
effectsize_bins
propde_bins

simList <- list()
for (effect in unique(MasterRes$effect)) {
  for (sampleType in unique(MasterRes$sampleType)) {
    for (m in unique(MasterRes$m)) {
      #for(sims in c("No Compensation","With Compensation")){
      for (sims in "No Compensation") {
        plotDflist <- list()
        plotlist <- list()
        ## Subset plots
        #subset <- MasterRes
        subset <-
          MasterRes[#!(as.character(MasterRes$test) %in% tests  & as.character(MasterRes$norm) %!in% norms) &
            #!(as.character(MasterRes$norm) == "doublerank") &
            #MasterRes$effect >=6 &
            #MasterRes$propde < .1 &
            MasterRes$test_norm %in% test_norms &
              MasterRes$m == m &
              MasterRes$effect == effect &
              MasterRes$sim == sims &
              MasterRes$sampleType == sampleType &
              MasterRes$propde > 0, ]#&
        if (nrow(subset) == 0) {
          next
        }
        subset[, fillvar][subset[, fillvar] == "t-test_rank-simple"] <-
          "t-test_rank"
        subset$test_norm[subset$test_norm == "t-test_rank-simple"] <-
          "t-test_rank"
        subset <-
          subset[subset$test_norm %!in% c("t-test_rank-corrected", "wilcox_rank-corrected"), ]
        subset$fdr <- as.numeric(as.character(subset$fdr))
        subset$sens <- as.numeric(as.character(subset$sens))
        subset$test_norm <-
          ifelse(
            subset$test_norm == "rankUnifTest_rank",
            "t-test_rank-corrected",
            subset$test_norm
          )
        subset$test_norm <-
          ifelse(
            subset$test_norm == "wilcoxUnifTest_rank",
            "wilcox_rank-corrected",
            subset$test_norm
          )
        # subset$test <- ifelse(subset$test_norm == "rankUnifTest_rank",
        #                       "t-test",
        #                       subset$test)
        # subset$test <- ifelse(subset$test_norm == "wilcoxUnifTest_rank",
        #                       "wilcox",
        #                       subset$test)
        # subset$norm <- ifelse(subset$test_norm == "rankUnifTest_rank",
        #                       "rank-corrected",
        #                       subset$norm)
        # subset$norm <- ifelse(subset$test_norm == "wilcoxUnifTest_rank",
        #                       "rank-corrected",
        #                       subset$norm)
        subset$test_norm <- ifelse(subset$test_norm == "DESeq2_ratio",
                                   "DESeq2_rle",
                                   subset$test_norm)
        # subset$test_norm <- ifelse(subset$test_norm == "t-test_ratio",
        #                            "t-test_rle",
        #                            subset$test_norm)
        # subset$test_norm <- ifelse(subset$test_norm == "wilcox_ratio",
        #                            "wilcox_rle",
        #                            subset$test_norm)
        # subset$test_norm <- ifelse(subset$test_norm == "wilcox_rank",
        #                            "wilcox_rank-simple",
        #                           subset$test_norm)
        subset$test_norm <- ifelse(subset$test_norm == "aldexT_aldex",
                                   "ALDEx2",
                                   subset$test_norm)
        # length <- length(unique(subset[, fillvar]))
        # custom_colors <- colorRamps::cyan2yellow(length)
        # custom_colors[regexpr("wilcox", sort(unique(subset$test_norm))) > 0] <-
        #   colorRamps::blue2green(length * 4)[(length + 2):(1 + length + sum(regexpr("wilcox", sort(
        #     unique(subset$test_norm)
        #   )) > 0))]
        # custom_colors[regexpr("t-test", sort(unique(subset$test_norm))) > 0] <-
        #   colorRamps::blue2red(length * 4)[1:sum(regexpr("t-test", sort(unique(
        #     subset$test_norm
        #   ))) > 0)]
        # custom_colors[regexpr("rank", sort(unique(subset$test_norm))) > 0] <-
        #   c("gold", "gold", "gold", "gold")[1:sum(regexpr("rank", sort(unique(
        #     subset$test_norm
        #   ))) > 0)]
        # names(custom_colors) <- sort(unique(subset$test_norm))
        # custom_colors[c(
        #   "DESeq2_rle",
        #   "edgeR_tmm",
        #   "metagenomeSeq_css",
        #   "Voom_tmm",
        #   "t-test_aldex",
        #   "wilcox_aldex"
        # )] <- colorRamps::green2red(20)[c(3, 5, 7, 9, 11, 13)]
        # 
        # plotTable <-
        #   as.data.frame(table(subset[, xvar], subset[, fillvar], subset[, facet]))
        # plotTable <- plotTable[plotTable[, "Freq"] != 0, ]
        # 
        # 
        # ## our differential abundance tests
        # order <-
        #   c(
        #     "DESeq2_rle",
        #     "edgeR_tmm",
        #     "metagenomeSeq_css",
        #     "Voom_tmm",
        #     "t-test_aldex",
        #     "t-test_none",
        #     "t-test_tss",
        #     "t-test_css",
        #     "t-test_tmm",
        #     "wilcox_none",
        #     "wilcox_tss",
        #     "wilcox_css",
        #     "wilcox_tmm",
        #     "t-test_rank-simple"
        #   )
        length <- length(unique(subset[,fillvar]))
        custom_colors <- colorRamps::cyan2yellow(length)
        custom_colors[regexpr("wilcox", sort(unique(subset$test_norm))) > 0] <-
          colorRamps::blue2green(length*4)[(length+2):(1+length+sum(regexpr("wilcox", sort(unique(subset$test_norm))) > 0))]
        custom_colors[regexpr("t-test", sort(unique(subset$test_norm))) > 0] <-
          colorRamps::blue2red(length*4)[1:sum(regexpr("t-test", sort(unique(subset$test_norm))) > 0)]
        custom_colors[regexpr("rank", sort(unique(subset$test_norm))) > 0] <-
          c("gold","gold","gold","gold")[1:sum(regexpr("rank", sort(unique(subset$test_norm))) > 0)]
        custom_colors[regexpr("perm", sort(unique(subset$test_norm))) > 0] <-
          c("yellowgreen","yellowgreen","yellowgreen","yellowgreen")[1:sum(regexpr("perm", sort(unique(subset$test_norm))) > 0)]
        names(custom_colors) <- sort(unique(subset$test_norm))
        custom_colors[c(
          "DESeq2",
          "edgeR_tmm",
          "metagenomeSeq_css",
          "Voom_tmm",
          "ALDEx2",
          "wilcox_aldex")] <- colorRamps::green2red(20)[c(3,5,7,9,11,13)]
        plotTable <- as.data.frame(table(subset[,xvar], subset[,fillvar], subset[,facet]))
        plotTable <- plotTable[plotTable[,"Freq"] != 0,]
        order <- 
          c("DESeq2",
            "edgeR_tmm",
            "metagenomeSeq_css",
            "Voom_tmm",
            "ALDEx2",
            "wilcox_aldex",
            "t-test_none",
            "t-test_tss",
            "t-test_css",
            "t-test_gmpr",
            "t-test_tmm",
            "t-test_wrench",
            "wilcox_none",
            "wilcox_tss",
            "wilcox_css",
            "wilcox_gmpr",
            "wilcox_tmm",
            "wilcox_wrench",
            "t-test_rank-simple",
            "perm-test_rank",
            "t-test_rank-corrected",
            "t-test_rank-mid-norm",
            "t-test_rank-mid",
            "wilcox_rank-corrected",
            "wilcox_rank-simple")
        ind <- sapply(order, function(x)
          grep(x, plotTable$Var1))
        
        custom_colors <- custom_colors[as.character(plotTable$Var1)]
        p <-
          ggplot(data = data.frame(x = seq(0, 1, .1), y = seq(0, 1, .1)), aes(x, y)) +
          xlim(0, 1) + ylim(0, 1) + xlab("1-FDR") + ylab("SENS")
        
        plotDf <- data.frame()
        gr <- grep("t-test_rank", plotTable$Var1)
        lr <- nrow(plotTable)
        plotTable[c(gr, lr), ] <- plotTable[c(lr, gr), ]
        custom_colors <- custom_colors[as.character(plotTable$Var1)]
        
        plotTable <- plotTable[plotTable$Var1 %in% incl_test_norms, ]
        custom_colors <-
          custom_colors[names(custom_colors) %in% incl_test_norms]
        for (i in 1:nrow(plotTable)) {
          tempxvar <- as.character(plotTable[i, 1])
          tempfillvar <- as.character(plotTable[i, 2])
          tempfacet <- as.character(plotTable[i, 3])
          
          temp <-
            subset[tolower(as.character(subset[, xvar])) == tolower(tempxvar) &
                     as.character(subset[, fillvar]) %in% tempfillvar  &
                     as.character(subset[, facet]) %in% tempfacet, ]
          mean_sens <-
            max(min(mean(
              as.numeric(as.character(temp[, "sens"])), na.rm = TRUE
            ), 1), 0)
          median_sens <-
            median(as.numeric(as.character(temp[, "sens"])), na.rm = TRUE)
          sd_sens <-
            sd(as.numeric(as.character(temp[, "sens"])), na.rm = TRUE)
          plus_sd <- max(mean_sens - sd_sens, 0)
          minus_sd <- min(1, mean_sens + sd_sens)
          firstquant <-
            quantile(as.numeric(as.character(temp[, "sens"])), .25, na.rm = TRUE)
          thirdquant <-
            quantile(as.numeric(as.character(temp[, "sens"])), .75, na.rm = TRUE)
          mean_fdr <-
            max(min(mean(
              as.numeric(as.character(temp[, "fdr"])), na.rm = TRUE
            ), 1), 0)
          median_fdr <-
            median(as.numeric(as.character(temp[, "fdr"])), na.rm = TRUE)
          sd_fdr <-
            sd(as.numeric(as.character(temp[, "fdr"])), na.rm = TRUE)
          plus_sd <- max(mean_fdr - sd_fdr, 0)
          minus_sd <- min(1, mean_fdr + sd_fdr)
          # mean distance for pairings of fdr & sensitivity vs. mean point
          # this will be our radius for circles
          msqe <- mean(sqrt((mean_sens - temp$sens) ^ 2 +
                              (mean_fdr - temp$fdr) ^ 2))
          temp$one_minus_fdr <- 1 - temp$fdr
          tempDf <- data.frame(
            xvar = tempxvar,
            fillvar = tempfillvar,
            facet = tempfacet,
            mean_sens = mean_sens,
            mean_fdr = mean_fdr,
            median_sens = median_sens,
            sd_sens = sd_sens,
            median_fdr = median_fdr,
            sd_fdr = sd_fdr,
            one_minus_fdr = 1 - mean_fdr,
            msqe = msqe,
            #ymin = minus_sd,
            #ymax = plus_sd,
            ymin = firstquant,
            ymax = thirdquant,
            stringsAsFactors = FALSE
          )
          plotDf <- rbind(plotDf, tempDf)
          # Starting the new plot
          area <-
            data.frame(x = temp$one_minus_fdr[chull(temp$one_minus_fdr,
                                                    temp$sens)],
                       y = temp$sens[chull(temp$one_minus_fdr,
                                           temp$sens)])
          point <-
            data.frame(x = tempDf$one_minus_fdr,
                       y = tempDf$mean_sens)
          alpha <- ifelse(tempDf$fillvar == "t-test_rank",
                          .8,
                          .2)
          color = custom_colors[as.character(plotTable[i, "Var1"])]
          r1 <- tempDf$sd_fdr
          r2 <- tempDf$sd_sens
          xc <- tempDf$one_minus_fdr
          yc <- tempDf$median_sens
          p <- p + geom_bubble(
            r1 = r1,
            r2 = r2,
            xc = xc,
            yc = yc,
            color = color,
            fill = color,
            alpha = alpha
          )
        }
        plotDf <- plotDf[order(plotDf$fillvar), ]
        custom_colors <- custom_colors[order(names(custom_colors))]
        text_colors <-  rep("grey1", length(unique(plotDf$fillvar)))
        text_colors[-c(1, 2, 3)] <- "grey90"
        text_colors[grep("t-test_rank", plotDf$fillvar)] <- "black"
        text_sizes <- rep(12, length(text_colors))
        text_sizes[grep("t-test_rank", plotDf$fillvar)] <- 14
        plotTable <- plotTable[order(plotTable$Var1), ]
        point <- data.frame(
          x = plotDf$one_minus_fdr,
          y = plotDf$mean_sens,
          fill = custom_colors,
          color = custom_colors,
          stringsAsFactors = FALSE
        )
        point$test_norm <- rownames(point)
        #p
        mean_fdr <- plotDf$mean_fdr[plotDf$fillvar == "t-test_rank"]
        mean_sens <- plotDf$mean_sens[plotDf$fillvar == "t-test_rank"]
        if(any(grepl("perm-test_rank",plotDf$fillvar))){
          mean_fdr_perm <- plotDf$mean_fdr[plotDf$fillvar == "perm-test_rank"]
          mean_sens_perm <- plotDf$mean_sens[plotDf$fillvar == "perm-test_rank"]
        }
        
        
        
        p <- p + geom_point(data = point,
                            size = 30,
                            aes(
                              x = x,
                              y = y,
                              color = test_norm
                            )) +
          theme(
            #axis.text.x=element_blank(),
            axis.ticks.x = element_blank(),
            text = element_text(size = 60),
            legend.text = element_text(size = 90),
            axis.text.x = element_text(size = 75),
            axis.text.y = element_text(size = 75),
            axis.title = element_text(size = 85, face = "bold"),
            plot.title = element_text(
              hjust = 0.5,
              size = 95,
              face = "bold"
            )
          ) +
          geom_text(
            data = point,
            label = 1:nrow(plotDf),
            size = (text_sizes * 2),
            #color = c(rep("grey90", 3), rep("grey10",4), "grey90", rep("grey10", length(unique(plotDf$fillvar))-10))) +
            color = text_colors
          ) +
          guides(color = guide_legend(title = "")) +
          scale_color_manual(values = custom_colors,
                             labels = paste0(names(custom_colors),
                                             " [", 1:length, "]"),
          ) +
          geom_hline(yintercept = 0,
                     colour = "black",
                     size = .4) +
          geom_vline(xintercept = 0,
                     colour = "black",
                     size = .4) +
          geom_vline(
            xintercept = 0.95,
            colour = "green",
            size = 4,
            linetype = "longdash"
          ) +
          geom_vline(
            xintercept = 1 - mean_fdr,
            linetype = "twodash",
            colour = "gold4",
            size = 5
          ) +
          geom_hline(
            yintercept = mean_sens,
            linetype = "twodash",
            colour = "gold4",
            size = 5
          ) +
          # geom_vline(
          #   xintercept = 1 - mean_fdr_perm,
          #   linetype = "twodash",
          #   colour = "yellowgreen",
          #   size = 5
          # ) +
          # geom_hline(
          #   yintercept = mean_sens_perm,
          #   linetype = "twodash",
          #   colour = "yellowgreen",
          #   size = 5
          # ) +
          # labs(y = "SENSITIVITY", x = "1-FDR", fill = fillvar,
          #      title = paste0("M = ", m)) +
          labs(
            y = ifelse(as.numeric(effect) %in% c(1.5),
                       paste0("M = ", m),
                       ""),
            x = "",
            fill = fillvar,
            title = ifelse(as.numeric(m) %in% c(5,50),
                           paste0("Fold-Change = ", effect),
                           "")
          ) +
          theme(
            panel.grid.minor = element_line(colour = "grey80", size = 0.4),
            panel.grid.major = element_line(colour = "grey60", size = .5)
          ) +
          scale_y_continuous(
            minor_breaks = seq(0 , 1, .05),
            breaks = seq(0, 1, .10),
            limits = c(0, 1)
          ) +
          scale_x_continuous(
            minor_breaks = seq(0 , 1, .05),
            breaks = seq(0, 1, .10),
            limits = c(0, 1)
          )
        simList[[paste0(sims, m, sampleType, effect)]] <- p
        
        #write.csv(plotDf, file = paste0(sims,m,sampleType,effect, "NegBin.csv"))
        write.csv(plotDf,
                  file = paste0(sims, 
                                m, 
                                sampleType, 
                                effect, 
                                distr_label,
                                ".csv"))
      }
    }
  }
}
plotlist2 <- simList

hmp_small <-
  names(plotlist2)[c(
    grep('No Compensation5HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation15HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation25HMP Mid Vagina1.5', names(plotlist2))
  )]
hmp_large <-
  names(plotlist2)[c(
    grep('No Compensation50HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation100HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation150HMP Mid Vagina1.5', names(plotlist2))
  )]



hmp_3 <-
  names(plotlist2)[c(
    grep('No Compensation5HMP Mid Vagina3', names(plotlist2)),
    grep('No Compensation15HMP Mid Vagina3', names(plotlist2)),
    grep('No Compensation25HMP Mid Vagina3', names(plotlist2)),
    grep('No Compensation50HMP Mid Vagina3', names(plotlist2)),
    grep('No Compensation100HMP Mid Vagina3', names(plotlist2)),
    grep('No Compensation150HMP Mid Vagina3', names(plotlist2))
  )]
hmp_5 <-
  names(plotlist2)[c(
    grep('No Compensation5HMP Mid Vagina5', names(plotlist2)),
    grep('No Compensation15HMP Mid Vagina5', names(plotlist2)),
    grep('No Compensation25HMP Mid Vagina5', names(plotlist2)),
    grep('No Compensation50HMP Mid Vagina5', names(plotlist2)),
    grep('No Compensation100HMP Mid Vagina5', names(plotlist2)),
    grep('No Compensation150HMP Mid Vagina5', names(plotlist2))
  )]
hmp_1.5 <-
  names(plotlist2)[c(
    grep('No Compensation5HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation15HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation25HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation50HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation100HMP Mid Vagina1.5', names(plotlist2)),
    grep('No Compensation150HMP Mid Vagina1.5', names(plotlist2))
  )]
## Prepare Plots
pl <- c(plotlist2[hmp_small], plotlist2[hmp_large])
pl <- c(pl[c(grep("n5H", names(pl)),
             grep("n50H", names(pl)),
             grep("n15H", names(pl)),
             grep("n100H", names(pl)),
             grep("n25H", names(pl)),
             grep("n150H", names(pl)))])
gg <- ggarrange(
  plotlist = pl,
  ncol = 2,
  nrow = 3,
  vjust = .75,
  common.legend = TRUE,
  legend = "bottom"
)
MeanPlot <- annotate_figure(
  gg,
  top = text_grob(
    paste0(
      ifelse(distr_label == "NegBin", 
             "Negative Binomial Distribution",
             "Beta-Binomial Distribution")),
    size = 100,
    face = "bold"
  ),
  left = text_grob(
    "Sensitivity",
    size = 100,
    rot = 90,
    face = "bold",
    hjust = 1,
    vjust = 1
  ),
  bottom = text_grob("1-FDR", size = 100, face = "bold")
)

## Plot it 

#x11()
MeanPlot
ggsave(
  MeanPlot,
  file = paste0("MainSimulationResults1.5_",distr_label,".pdf"),
  width = 75,
  height = 79,
  limitsize = FALSE,
  device = "pdf"
)
# width = 75, # 75, 79 before
# height = 79,
# limitsize = FALSE,
# dpi = 600,
# device = "pdf")
#}


## Create Tables For Simulation Results

try({
df1 <- read.csv(paste0(datawd,"No Compensation50HMP Mid Vagina3",distr_label,".csv"))
df2 <- read.csv(paste0(datawd,"No Compensation50HMP Mid Vagina5",distr_label,".csv"))
df3 <- read.csv(paste0(datawd,"No Compensation100HMP Mid Vagina3",distr_label,".csv"))
df4 <- read.csv(paste0(datawd,"No Compensation100HMP Mid Vagina5",distr_label,".csv"))
df5 <- read.csv(paste0(datawd,"No Compensation150HMP Mid Vagina3",distr_label,".csv"))
df6 <- read.csv(paste0(datawd,"No Compensation150HMP Mid Vagina5",distr_label,".csv"))},silent = TRUE)

try({
df1 <- read.csv(paste0(datawd,"No Compensation5HMP Mid Vagina3",distr_label,".csv"))
df2 <- read.csv(paste0(datawd,"No Compensation5HMP Mid Vagina5",distr_label,".csv"))
df3 <- read.csv(paste0(datawd,"No Compensation15HMP Mid Vagina3",distr_label,".csv"))
df4 <- read.csv(paste0(datawd,"No Compensation15HMP Mid Vagina5",distr_label,".csv"))
df5 <- read.csv(paste0(datawd,"No Compensation25HMP Mid Vagina3",distr_label,".csv"))
df6 <- read.csv(paste0(datawd,"No Compensation25HMP Mid Vagina5",distr_label,".csv"))},silent = TRUE)

try({
  df1 <- read.csv(paste0(datawd,"No Compensation5HMP Mid Vagina1.5",distr_label,".csv"))
  df2 <- read.csv(paste0(datawd,"No Compensation15HMP Mid Vagina1.5",distr_label,".csv"))
  df3 <- read.csv(paste0(datawd,"No Compensation25HMP Mid Vagina1.5",distr_label,".csv"))
  df4 <- read.csv(paste0(datawd,"No Compensation50HMP Mid Vagina1.5",distr_label,".csv"))
  df5 <- read.csv(paste0(datawd,"No Compensation100HMP Mid Vagina1.5",distr_label,".csv"))
  df6 <- read.csv(paste0(datawd,"No Compensation150HMP Mid Vagina1.5",distr_label,".csv"))},silent = TRUE)

clean <- function(df, m, efx) {
  d <- data.frame(
    "Sensitivity" = paste0(round(df$mean_sens, 3),
                           " +/- (",
                           round(df$sd_sens, 3),
                           ")"),
    "FDR" = paste0(round(df$mean_fdr, 3),
                   " +/- (",
                   round(df$sd_fdr, 3),
                   ")"),
    "M" = m,
    "FC" = efx
    
  )
  rownames(d) <- as.character(df$fillvar)
  d
}
min <- min(MasterRes$m)
med <- mean(MasterRes$m)
max <- max(MasterRes$m)
low <- min(MasterRes$effect)
high <- max(MasterRes$effect)
df1 <- clean(df1, 5, low)
df2 <- clean(df2, 15, high)
df3 <- clean(df3, 25, low)
df4 <- clean(df4, 50, high)
df5 <- clean(df5, 100, low)
df6 <- clean(df6, 150, high)
df <- cbind(df1, df2, df3, df4, df5, df6)
rownames(df) <- sapply(rownames(df), function(x) {
  s <- strsplit(x, "_")[[1]]
  s1 <- s[1]
  s2 <- s[2]
  s1 <- paste0(toupper(substr(s1, 1, 1)),
               substr(s1, 2, nchar(s1)))
  s2 <- toupper(s2)
  paste0(s2, " ", s1)
})
rownames(df)[c("ALDEX")]
rownames(df)[grep("t-test_rank", rownames(df))] <- "Ranked T-Test"
rownames(df)[grep("RANK Perm-test", rownames(df))] <- "RANKED Permutation T-test"
df <- apply(df, 2, function(x) {
  ifelse(x == "0 +/- (0)",
         "0",
         x)
})
n <- ncol(df) / 4
df2 <- data.frame()
for (cols in 1:n) {
  df2 <- rbind(df2, cbind(rownames(df), df[, 1:4 + (cols - 1) * 4]))
}
colnames(df2)[1] <- "DAT"
colnames(df2) <-
  c(
    "Methodology",
    "Mean Sensitivity (+/- SD)",
    "Mean FDR (+/- SD)",
    "Sample Size Per Condition",
    "Fold-Change"
  )
#df2 <- df2[,1:5]
write.csv(df2, file = paste0("PlotTables21.5_",distr_label,".csv"), row.names = F)