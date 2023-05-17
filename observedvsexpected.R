dirs <- list.dirs("D:/exeter_NT001/", recursive = FALSE)
dirs <- dirs[!grepl("RECYCLE", dirs)]
print(dirs)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
for (d in dirs) {
  setwd("D:/H7")
  list29903 <- read_table("29903.csv")
  all <- read_table("underscore.txt")
  #print(all)
  selectedRows <- all[grep(substr(d, 18, 20), all$plate_well), ]
  #print(selectedRows)
  alphaEF <- as.numeric(selectedRows$alpha_exp_frequency[1])
  #print(alphaEF)
  betaEF <- selectedRows$beta_exp_frequency[1]
  deltaC23EF <- selectedRows$deltaC23_exp_frequency[1]
  deltaC29EF <- selectedRows$deltaC29_exp_frequency[1]
  omicronEF <- selectedRows$omicron_exp_frequency[1]
  alphacontrol <- read_table("AlphaC15.tsv")
  alphacontrol$expectedfrequency <- alphaEF
  alphacontrol$lineage <- c("alpha")
  betacontrol <- read_table("BetaC16.tsv")
  betacontrol$expectedfrequency <- betaEF
  betacontrol$lineage <- c("beta")
  deltaC23control <- read_table("DeltaC23.tsv")
  deltaC23control$expectedfrequency <- deltaC23EF
  deltaC23control$lineage <- c("deltaC23")
  deltaC29control <- read_table("Delta-AY.2-C29.tsv")
  deltaC29control$expectedfrequency <- deltaC29EF
  deltaC29control$lineage <- c("deltaC29")
  omicroncontrol <- read_table("Omicron-C48.tsv")
  omicroncontrol$expectedfrequency <- omicronEF
  omicroncontrol$lineage <- c("omicron")
  expfrequencyonly <- select(selectedRows, c(alpha_exp_frequency,beta_exp_frequency, deltaC23_exp_frequency, deltaC29_exp_frequency, omicron_exp_frequency))
  #print(expfrequencyonly)
  expfrequencyonly <- expfrequencyonly[, colSums(expfrequencyonly != 0) > 0]
  list <- colnames(expfrequencyonly)
  list <- gsub('_exp_frequency','', list)
  List2 <- paste(unlist(list),collapse="")
  #print(List2)
  # Set the working directory
  setwd(d)
  # Read the LR output file
  files <- list.files(pattern = "\\.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv$")
  #print(files)
  LRoutput <- read.table(files[1], header=TRUE)
  LRoutput <- separate(data = LRoutput, col = "Cons.Cov.Reads1.Reads2.Freq.P.value", into = c("Col", "Cov", "Reads1", "Reads2", "VarFreq", "P-value"), sep = "\\:")
  LRoutput$VarFreq <- as.numeric(sub("%", "", LRoutput$VarFreq, fixed = TRUE)) / 100
  #print(LRoutput)
  # Read the LR depths file
  depth <- list.files(pattern = "\\.ont.samtools.cov")
  LRdepths <- read.table(depth[1], col.names = c("NAME", "POS", "DEPTH"))
  allcontrols <- rbind(alphacontrol, betacontrol, deltaC23control, deltaC29control, omicroncontrol)
  total <- aggregate(cbind(expectedfrequency) ~ POS, data = allcontrols, sum)
  truecontrols <- rbind(alphacontrol, betacontrol, deltaC23control, deltaC29control, omicroncontrol)
  library(dplyr)
  realSNPsforfalse <- truecontrols[truecontrols$lineage %in% list,]
  LRvsTRUE <- merge(LRoutput, realSNPsforfalse, by.x=c("Position", "Var"), by.y=c("POS", "ALT"))
  TP <- as.data.frame(LRvsTRUE)
  FP <- LRoutput %>% anti_join(realSNPsforfalse, by = c("Position" = "POS"))
  FN <- realSNPsforfalse %>% anti_join(LRvsTRUE, by = c("POS" = "Position"))
  LRoutput$DEPTH <- LRoutput$Cov
  FN$DEPTH <- LRdepths$DEPTH[match(FN$POS, LRdepths$POS)]
  FNbecausecoverage <- FN[FN$DEPTH < 3,]
  FNbecausejustnotcalled <- FN[FN$DEPTH > 3,]
  TNBEFORELR <- list29903 %>% anti_join(realSNPsforfalse, by = c("POS" ="POS"))
  TN <- TNBEFORELR %>% anti_join(LRoutput, by = c("POS" = "Position"))
  total$observedfrequency <- LRoutput$VarFreq[match(total$POS, LRoutput$Position)]
  total$depth <- LRdepths$DEPTH[match(total$POS, LRdepths$POS)]
  total$lineage <- character(nrow(total))
  deltaC23_freq <- unlist(selectedRows$`deltaC23_exp_frequency`)
  total$deltaC23_freq <- deltaC23_freq[1]
  omicron_conc <- unlist(selectedRows$`omicron-C48`)
  total$omicron_conc <- omicron_conc[1]
  for (i in 1:nrow(total)) {
    matching_indices <- which(realSNPsforfalse$POS == total$POS[i])
    if (length(matching_indices) == 0) {
      total$lineage[i] <- "Control SNPs"
    } else if (length(matching_indices) == 1) {
      total$lineage[i] <- realSNPsforfalse$lineage[matching_indices]
    } else {
      total$lineage[i] <- paste(realSNPsforfalse$lineage[matching_indices], collapse = ",")
    }
  }
  total[is.na(total)] <- 0
  total$sample <- substr(d,18,20)
  assign(paste0("total_ONT", substr(d, 18,20)), total)
  # Define the list of controls
}



all3ONT <- do.call(rbind, lapply( ls(patt="total_ONT"), get) )
all3ONT$sequencing <- "ONT"
all3ONT <- all3ONT[all3ONT$sample == "A1_" |all3ONT$sample == "B1_" |all3ONT$sample == "C1_" | all3ONT$sample == "D1_" |all3ONT$sample == "E1_",]
all3ONT$color <- "#FF61CC"
all3ONT$color <- ifelse(grepl("alpha$", all3ONT$lineage), "#F8766D", all3ONT$color)
all3ONT$color <- ifelse(grepl("alpha,beta$", all3ONT$lineage), "#762A01", all3ONT$color)
all3ONT$color <- ifelse(grepl("^beta", all3ONT$lineage), "#CD9600", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC23$", all3ONT$lineage), "#7CAE00", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC29$", all3ONT$lineage), "#0E5615", all3ONT$color)
all3ONT$color <- ifelse(grepl("omicron$", all3ONT$lineage), "#A58AFF", all3ONT$color)
all3ONT$color <- as.factor(all3ONT$color)
all3ONT <- all3ONT %>% 
  group_by(sequencing) %>% 
  mutate(normalized_depth = depth / max(depth))
p <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, color=color, alpha=normalized_depth))
g <- p + geom_point(stat="identity", color=all3ONT$color, alpha=all3ONT$normalized_depth, size = 2.5) + xlab("Expected Frequency") + ylab("Observed Frequency")
final9 <- g+ scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#762A01" = "alpha,beta", "#CD9600" = "beta", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.3,1.0))
corr <- cor.test(x=all3ONT$expectedfrequency, y=all3ONT$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr
ggsave("OBSvsEXP_EXETER_NT001_CovMix1-5_ONT.png", final9)

my_hist <- ggplot(all3ONT, aes(expectedfrequency, fill = all3ONT$color)) + 
  geom_bar() 
my_hist
snplegend <- my_hist + scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#762A01" = "alpha,beta", "#CD9600" = "beta", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.1,1.0))
snplegend
library(grid)
legend <- cowplot::get_legend(snplegend)
grid.newpage()
grid.draw(legend)
ggsave("CovMix1-5_ONT_SNPs_legend.png", legend)

my_hist2 <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, alpha = all3ONT$normalized_depth)) + 
  geom_point(size=3)
my_hist2
snplegend2 <- my_hist2 + scale_alpha_continuous(name = "Normalized ONT Depth")
snplegend2


library(grid)
legend2 <- cowplot::get_legend(snplegend2)
grid.newpage()
grid.draw(legend2)
ggsave("OBSvsEXP_EXETER_NT001_CovMix1-5_ONT_alpha_legend.png",legend2)





all3ONT <- do.call(rbind, lapply( ls(patt="total_ONT"), get) )
all3ONT$sequencing <- "ONT"
all3ONT <- all3ONT[all3ONT$sample == "F1_" |all3ONT$sample == "G1_" |all3ONT$sample == "H1_" | all3ONT$sample == "A2_" |all3ONT$sample == "B2_" |all3ONT$sample == "C2_" |all3ONT$sample == "D2_" |all3ONT$sample == "E2_" |all3ONT$sample == "F2_" |all3ONT$sample == "G2_" |all3ONT$sample == "H2_"     ,]
all3ONT$color <- "#FF61CC"
all3ONT$color <- ifelse(grepl("alpha$", all3ONT$lineage), "#F8766D", all3ONT$color)
all3ONT$color <- ifelse(grepl("alpha,beta$", all3ONT$lineage), "#762A01", all3ONT$color)
all3ONT$color <- ifelse(grepl("alpha,beta,deltaC23$", all3ONT$lineage), "#353302", all3ONT$color)

all3ONT$color <- ifelse(grepl("^beta", all3ONT$lineage), "#CD9600", all3ONT$color)
all3ONT$color <- ifelse(grepl("beta,deltaC23$", all3ONT$lineage), "#4A7601", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC23$", all3ONT$lineage), "#7CAE00", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC29$", all3ONT$lineage), "#0E5615", all3ONT$color)
all3ONT$color <- ifelse(grepl("omicron$", all3ONT$lineage), "#A58AFF", all3ONT$color)
all3ONT$color <- as.factor(all3ONT$color)
all3ONT <- all3ONT %>% 
  group_by(sequencing) %>% 
  mutate(normalized_depth = depth / max(depth))
p <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, color=color, alpha=normalized_depth))
g <- p + geom_point(stat="identity", color=all3ONT$color, alpha=all3ONT$normalized_depth, size = 2.5) + xlab("Expected Frequency") + ylab("Observed Frequency")
final9 <- g+ scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#762A01" = "alpha,beta", "#7CAE00" = "deltaC23", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.3,1.0))

corr <- cor.test(x=all3ONT$expectedfrequency, y=all3ONT$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr
ggsave("OBSvsEXP_EXETER_NT001_CovMix6-16_ONT.png", final9)

my_hist <- ggplot(all3ONT, aes(expectedfrequency, fill = all3ONT$color)) + 
  geom_histogram() 
my_hist
snplegend <- my_hist + scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#762A01" = "alpha,beta", "#CD9600" = "beta", "#4A7601" = "beta,deltaC23", "#353302" = "alpha,beta,deltaC23", "#7CAE00" = "deltaC23", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.1,1.0))
snplegend
library(grid)
legend <- cowplot::get_legend(snplegend)
grid.newpage()
grid.draw(legend)
ggsave("CovMix6-16_ONT_SNPs_legend.png", legend)

my_hist2 <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, alpha = all3ONT$normalized_depth)) + 
  geom_point(size=3)
my_hist2
snplegend2 <- my_hist2 + scale_alpha_continuous(name = "Normalized ONT Depth")
snplegend2


library(grid)
legend2 <- cowplot::get_legend(snplegend2)
grid.newpage()
grid.draw(legend2)
ggsave("OBSvsEXP_EXETER_NT001_CovMix6-16_ONT_alpha_legend.png",legend2)

all3ONT <- do.call(rbind, lapply( ls(patt="total_ONT"), get) )
all3ONT$sequencing <- "ONT"
all3ONT <- all3ONT[all3ONT$sample == "A3_" |all3ONT$sample == "B3_" |all3ONT$sample == "C3_" | all3ONT$sample == "D3_" |all3ONT$sample == "E3_" |all3ONT$sample == "F3_" |all3ONT$sample == "G3_" |all3ONT$sample == "A4_" |all3ONT$sample == "B4_" |all3ONT$sample == "C4_" |all3ONT$sample== "D4_" |all3ONT$sample == "E4_" |all3ONT$sample == "F4_"  ,]
all3ONT$color <- "#FF61CC"
all3ONT$color <- ifelse(grepl("alpha$", all3ONT$lineage), "#F8766D", all3ONT$color)
all3ONT$color <- ifelse(grepl("^beta", all3ONT$lineage), "#CD9600", all3ONT$color)
all3ONT$color <- ifelse(grepl("beta,deltaC23$", all3ONT$lineage), "#4A7601", all3ONT$color)
all3ONT$color <- ifelse(grepl("beta,deltaC29$", all3ONT$lineage), "#0E2702", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC23$", all3ONT$lineage), "#7CAE00", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC29$", all3ONT$lineage), "#0E5615", all3ONT$color)
all3ONT$color <- ifelse(grepl("^deltaC23,deltaC29", all3ONT$lineage), "#063970", all3ONT$color)
all3ONT$color <- ifelse(grepl("^beta,deltaC23,deltaC29$", all3ONT$lineage), "#050F01", all3ONT$color)
all3ONT$color <- ifelse(grepl("omicron$", all3ONT$lineage), "#A58AFF", all3ONT$color)
all3ONT$color <- as.factor(all3ONT$color)
all3ONT <- all3ONT %>% 
  group_by(sequencing) %>% 
  mutate(normalized_depth = depth / max(depth))
p <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, color=color, alpha=normalized_depth))
g <- p + geom_point(stat="identity", color=all3ONT$color, alpha=all3ONT$normalized_depth, size = 2.5) + xlab("Expected Frequency") + ylab("Observed Frequency")
final9 <- g+ scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#0E5615" = "deltaC29", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.3,1.0))
corr <- cor.test(x=all3ONT$expectedfrequency, y=all3ONT$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr
ggsave("OBSvsEXP_EXETER_NT001_CovMix17-29_ONT.png", final9)

my_hist <- ggplot(all3ONT, aes(expectedfrequency, fill = all3ONT$color)) + 
  geom_histogram() 
my_hist
snplegend <- my_hist + scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#050F01" = "beta,deltaC23,deltaC29", "#7CAE00" = "deltaC23", "#063970" = "deltaC23,deltaC29", "#0E5615" = "deltaC29", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.1,1.0))
snplegend
library(grid)
legend <- cowplot::get_legend(snplegend)
grid.newpage()
grid.draw(legend)
ggsave("CovMix17-29_ONT_SNPs_legend.png", legend)

my_hist2 <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, alpha = all3ONT$normalized_depth)) + 
  geom_point(size=3)
my_hist2
snplegend2 <- my_hist2 + scale_alpha_continuous(name = "Normalized ONT Depth")
snplegend2


library(grid)
legend2 <- cowplot::get_legend(snplegend2)
grid.newpage()
grid.draw(legend2)
ggsave("OBSvsEXP_EXETER_NT001_CovMix17-29_ONT_alpha_legend.png",legend2)

all3ONT <- do.call(rbind, lapply( ls(patt="total_ONT"), get) )
all3ONT$sequencing <- "ONT"
all3ONT <- all3ONT[all3ONT$sample == "B7_" |all3ONT$sample== "D7_" |all3ONT$sample == "E7_" |all3ONT$sample == "F7_" | all3ONT$sample == "G7_" | all3ONT$sample == "H7_" | all3ONT$sample == "C7_" ,]
all3ONT$color <- "#FF61CC"
all3ONT$color <- ifelse(grepl("alpha$", all3ONT$lineage), "#F8766D", all3ONT$color)
all3ONT$color <- ifelse(grepl("beta$", all3ONT$lineage), "#CD9600", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC23$", all3ONT$lineage), "#7CAE00", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC29$", all3ONT$lineage), "#0E5615", all3ONT$color)
all3ONT$color <- ifelse(grepl("omicron$", all3ONT$lineage), "#A58AFF", all3ONT$color)
all3ONT$color <- as.factor(all3ONT$color)
all3ONT <- all3ONT %>% 
  group_by(sequencing) %>% 
  mutate(normalized_depth = depth / max(depth))
p <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, color=color, alpha=normalized_depth))
g <- p + geom_point(stat="identity", color=all3ONT$color, alpha=all3ONT$normalized_depth, size = 2.5) + xlab("Expected Frequency") + ylab("Observed Frequency")
final9 <- g+ scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#0E5615" = "deltaC29", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.3,1.0))
corr <- cor.test(x=all3ONT$expectedfrequency, y=all3ONT$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr
ggsave("OBSvsEXP_EXETER_NT001_Alpha_sd_1-7_ONT.png", final9)

my_hist <- ggplot(all3ONT, aes(expectedfrequency, fill = all3ONT$color)) + 
  geom_histogram() 
my_hist
snplegend <- my_hist + scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#0E5615" = "deltaC29", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.1,1.0))
snplegend
library(grid)
legend <- cowplot::get_legend(snplegend)
grid.newpage()
grid.draw(legend)
ggsave("Alpha_SD1-7_ONT_SNPs_legend.png", legend)

my_hist2 <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, alpha = all3ONT$normalized_depth)) + 
  geom_point(size=3)
my_hist2
snplegend2 <- my_hist2 + scale_alpha_continuous(name = "Normalized ONT Depth")
snplegend2


library(grid)
legend2 <- cowplot::get_legend(snplegend2)
grid.newpage()
grid.draw(legend2)
ggsave("OBSvsEXP_EXETER_NT001_Alpha_sd_1-7_ONT_alpha_legend.png",legend2)

all3ONT <- do.call(rbind, lapply( ls(patt="total_ONT"), get) )
all3ONT$sequencing <- "ONT"
all3ONT <- all3ONT[all3ONT$sample == "B8_" |all3ONT$sample== "D8_" |all3ONT$sample == "E8_" |all3ONT$sample == "F8_" | all3ONT$sample == "G8_" | all3ONT$sample == "H8_" | all3ONT$sample == "C8_" ,]
all3ONT$color <- "#FF61CC"
all3ONT$color <- ifelse(grepl("alpha$", all3ONT$lineage), "#F8766D", all3ONT$color)
all3ONT$color <- ifelse(grepl("beta$", all3ONT$lineage), "#CD9600", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC23$", all3ONT$lineage), "#7CAE00", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC29$", all3ONT$lineage), "#0E5615", all3ONT$color)
all3ONT$color <- ifelse(grepl("omicron$", all3ONT$lineage), "#A58AFF", all3ONT$color)
all3ONT$color <- as.factor(all3ONT$color)
all3ONT <- all3ONT %>% 
  group_by(sequencing) %>% 
  mutate(normalized_depth = depth / max(depth))
p <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, color=color, alpha=normalized_depth))
g <- p + geom_point(stat="identity", color=all3ONT$color, alpha=all3ONT$normalized_depth, size = 2.5) + xlab("Expected Frequency") + ylab("Observed Frequency")
final9 <- g+ scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#0E5615" = "deltaC29", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.3,1.0))
final9
corr <- cor.test(x=all3ONT$expectedfrequency, y=all3ONT$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr
ggsave("OBSvsEXP_EXETER_NT001_Delta_sd_1-7_ONT.png", final9)

my_hist <- ggplot(all3ONT, aes(expectedfrequency, fill = all3ONT$color)) + 
  geom_histogram() 
my_hist
snplegend <- my_hist + scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#0E5615" = "deltaC29", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.1,1.0))
snplegend
library(grid)
legend <- cowplot::get_legend(snplegend)
grid.newpage()
grid.draw(legend)
ggsave("Delta_SD1-7_ONT_SNPs_legend.png", legend)

my_hist2 <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, alpha = all3ONT$normalized_depth)) + 
  geom_point(size=3)
my_hist2
snplegend2 <- my_hist2 + scale_alpha_continuous(name = "Normalized ONT Depth")
snplegend2


library(grid)
legend2 <- cowplot::get_legend(snplegend2)
grid.newpage()
grid.draw(legend2)
ggsave("OBSvsEXP_EXETER_NT001_Delta_sd_1-7_ONT_alpha_legend.png",legend2)

all3ONT <- do.call(rbind, lapply( ls(patt="total_ONT"), get) )
all3ONT$sequencing <- "ONT"
all3ONT <- all3ONT[all3ONT$sample == "B9_" |all3ONT$sample== "D9_" |all3ONT$sample == "E9_" |all3ONT$sample == "F9_" | all3ONT$sample == "G9_" | all3ONT$sample == "H9_" | all3ONT$sample == "C9_" ,]
all3ONT$color <- "#FF61CC"
all3ONT$color <- ifelse(grepl("alpha$", all3ONT$lineage), "#F8766D", all3ONT$color)
all3ONT$color <- ifelse(grepl("alpha,deltaC23$", all3ONT$lineage), "#3E1006", all3ONT$color)
all3ONT$color <- ifelse(grepl("beta$", all3ONT$lineage), "#CD9600", all3ONT$color)
all3ONT$color <- ifelse(grepl("^deltaC23", all3ONT$lineage), "#7CAE00", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC29$", all3ONT$lineage), "#0E5615", all3ONT$color)
all3ONT$color <- ifelse(grepl("omicron$", all3ONT$lineage), "#A58AFF", all3ONT$color)
all3ONT$color <- as.factor(all3ONT$color)
all3ONT <- all3ONT %>% 
  group_by(sequencing) %>% 
  mutate(normalized_depth = depth / max(depth))
p <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, color=color, alpha=normalized_depth))
g <- p + geom_point(stat="identity", color=all3ONT$color, alpha=all3ONT$normalized_depth, size = 2.5) + xlab("Expected Frequency") + ylab("Observed Frequency")
final9 <- g+ scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#0E5615" = "deltaC29", "#9BA20B" = "alpha,deltaC23", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.3,1.0))
final9
corr <- cor.test(x=all3ONT$expectedfrequency, y=all3ONT$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr
ggsave("OBSvsEXP_EXETER_NT001_Alpha_Delta_sd_1-7_ONT.png", final9)

my_hist <- ggplot(all3ONT, aes(expectedfrequency, fill = all3ONT$color)) + 
  geom_histogram() 
my_hist
snplegend <- my_hist + scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#3E1006" = "alpha,deltaC23", "#0E5615" = "deltaC29", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.1,1.0))
snplegend
library(grid)
legend <- cowplot::get_legend(snplegend)
grid.newpage()
grid.draw(legend)
ggsave("Alpha_Delta_SD1-7_ONT_SNPs_legend.png", legend)

my_hist2 <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, alpha = all3ONT$normalized_depth)) + 
  geom_point(size=3)
my_hist2
snplegend2 <- my_hist2 + scale_alpha_continuous(name = "Normalized ONT Depth")
snplegend2


library(grid)
legend2 <- cowplot::get_legend(snplegend2)
grid.newpage()
grid.draw(legend2)
ggsave("OBSvsEXP_EXETER_NT001_Alpha_Delta_sd_1-7_ONT_alpha_legend.png",legend2)

all3ONT <- do.call(rbind, lapply( ls(patt="total_ONT"), get) )
all3ONT$sequencing <- "ONT"
all3ONT <- all3ONT[all3ONT$sample == "A11" |all3ONT$sample == "B11" |all3ONT$sample== "D11" |all3ONT$sample == "E11" |all3ONT$sample == "F11" | all3ONT$sample == "G11" | all3ONT$sample == "H11" | all3ONT$sample == "C11"| all3ONT$sample == "A12" |all3ONT$sample == "B12" | all3ONT$sample == "C12" | all3ONT$sample == "D12" | all3ONT$sample == "E12" ,]
all3ONT$color <- "#FF61CC"
all3ONT$color <- ifelse(grepl("alpha$", all3ONT$lineage), "#F8766D", all3ONT$color)
all3ONT$color <- ifelse(grepl("beta$", all3ONT$lineage), "#CD9600", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC23$", all3ONT$lineage), "#7CAE00", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC29$", all3ONT$lineage), "#0E5615", all3ONT$color)
all3ONT$color <- ifelse(grepl("deltaC23,omicron$", all3ONT$lineage), "#170854", all3ONT$color)
all3ONT$color <- ifelse(grepl("^omicron", all3ONT$lineage), "#A58AFF", all3ONT$color)
all3ONT$color <- as.factor(all3ONT$color)
all3ONT <- all3ONT %>% 
  group_by(sequencing) %>% 
  mutate(normalized_depth = depth / max(depth))
p <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, color=color, alpha=normalized_depth))
g <- p + geom_point(stat="identity", color=all3ONT$color, alpha=all3ONT$normalized_depth, size = 2.5) + xlab("Expected Frequency") + ylab("Observed Frequency")
final9 <- g+ scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#0E5615" = "deltaC29", "#170854" = "deltaC23,omicron", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.3,1.0))
final9
corr <- cor.test(x=all3ONT$expectedfrequency, y=all3ONT$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr
ggsave("OBSvsEXP_EXETER_NT001_Delta_Omicron_sd_1-13_ONT.png", final9, width = 5, height = 5)

my_hist <- ggplot(all3ONT, aes(expectedfrequency, fill = all3ONT$color)) + 
  geom_histogram() 
my_hist
snplegend <- my_hist + scale_fill_identity(name="colour", guide = 'legend', labels=c("#F8766D" = "alpha", "#CD9600" = "beta", "#7CAE00" = "deltaC23", "#0E5615" = "deltaC29", "#170854" = "deltaC23,omicron", "#A58AFF" = "omicron", "#FF61CC" = "Control SNPs")) + geom_abline(linetype = "dotted") + theme_classic() + scale_alpha(range=c(0.1,1.0))
snplegend
library(grid)
legend <- cowplot::get_legend(snplegend)
grid.newpage()
grid.draw(legend)
ggsave("Delta_Omicron_SD1-13_ONT_SNPs_legend.png", legend)

my_hist2 <- ggplot(all3ONT, aes(x=expectedfrequency, y=observedfrequency, alpha = all3ONT$normalized_depth)) + 
  geom_point(size=3)
my_hist2
snplegend2 <- my_hist2 + scale_alpha_continuous(name = "Normalized ONT Depth")
snplegend2


library(grid)
legend2 <- cowplot::get_legend(snplegend2)
grid.newpage()
grid.draw(legend2)
ggsave("OBSvsEXP_EXETER_NT001_Delta_Omicron_sd_1-13_ONT_alpha_legend.png",legend2)






