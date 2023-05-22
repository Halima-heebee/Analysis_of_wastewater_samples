dirs <- list.dirs("/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002/", recursive = FALSE)
dirs <- dirs[!grepl("RECYCLE", dirs)]
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
for (d in dirs) {
  setwd("/home/mbxha18/allcombinations/New_illumina/exeter/DAVE")
  list29903 <- read_table("29903.csv")
  all <- read_table("Copy of TWIST_synthetic_mixes_sample_sheet_150222.txt")
  #print(all)
  selectedRows <- all[grep(substr(d, 83, 85), all$plate_well), ]
  print(selectedRows)
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
  deltaC29control$lineage <- c("deltaAY2")
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
  LRoutput <- read.table(files[1], header = TRUE)
  LRoutput <- separate(data = LRoutput, col = "Cons.Cov.Reads1.Reads2.Freq.P.value", into = c("Col", "Cov", "Reads1", "Reads2", "VarFreq", "P-value"), sep = "\\:")
  LRoutput$VarFreq <- as.numeric(sub("%", "", LRoutput$VarFreq, fixed = TRUE)) / 100
  #print(LRoutput)
  # Read the LR depths file
  depth <- list.files(pattern = "\\.ont.samtools.cov")
  LRdepths <- read.table(depth[1], col.names = c("NAME", "POS", "DEPTH"))
  
  # Define the list of controls
  control_list <- list(alphacontrol, betacontrol, deltaC23control, deltaC29control, omicroncontrol)
  allcontrols <- do.call(rbind, control_list)
  #print(allcontrols)
  # Compute the total expected frequency
  total_expected_frequency <- aggregate(cbind(expectedfrequency) ~ POS, data = allcontrols, sum)
  #print(total_expected_frequency)
  substrings <- c("alpha", "beta", "deltaC23", "deltaC29", "omicron")
  matches <- str_extract_all(List2, paste(substrings, collapse = "|"))[[1]]
  list3 <- unique(matches)
  # Filter the controls for the true SNPs
  truecontrols <- allcontrols[allcontrols$lineage %in% list3,]
  #print(truecontrols)  
  # Compute the TP and FP for the LR output
  LRvsTRUE <- merge(LRoutput, truecontrols, by.x = c("Position", "Var"), by.y = c("POS", "ALT"))
  TP_ONT <- as.data.frame(LRvsTRUE)
  assign(paste0("TP_ONT", substr(d, 11,100)), TP_ONT)
  
  
  FP_ONT <- anti_join(LRoutput, truecontrols, by = c("Position" = "POS"))
  assign(paste0("FP_ONT", substr(d, 11,100)), FP_ONT)
  # Compute the FN for the true SNPs
  FN_ONT <- anti_join(truecontrols, LRvsTRUE, by = c("POS" = "Position"))
  
  FN_ONT$DEPTH <- LRdepths$DEPTH[match(FN_ONT$POS, LRdepths$POS)]
  assign(paste0("FN_ONT", substr(d, 11,100)), FN_ONT)
  # Filter the FN based on coverage
  #FNbecausecoverageONT <- FN_ONT[[d]][FN_ONT[[d]]$DEPTH < 3,]
  #FNbecausejustnotcalledONT <- FN_ONT[[d]][FN_ONT[[d]]$DEPTH > 3,]
  
  # Compute the TN for the LR output
  TNBEFORELR <- anti_join(list29903, truecontrols, by = c("POS" = "POS"))
  TN_ONT <- anti_join(TNBEFORELR, LRoutput, by = c("POS" = "Position"))
  assign(paste0("TN_ONT", substr(d, 11,100)), TN_ONT)
  # Add the observed frequency, depth, lineage, and sample columns to the total expected frequency data frame
  total_expected_frequency$observedfrequency <- LRoutput$VarFreq[match(total_expected_frequency$POS, LRoutput$Position)]
  total_expected_frequency$depth <- LRdepths$DEPTH[match(total_expected_frequency$POS, LRdepths$POS)]
  total_expected_frequency$lineage <- List2
  total_expected_frequency[is.na(total_expected_frequency)] <- 0
  total_expected_frequency$sample <- d
  
  # Append the total expected frequency data frame to the list
  total_ONT <- total_expected_frequency
  assign(paste0("total_ONT", substr(d, 11,100)), total_ONT)
}

for (d in dirs) {
  setwd("/home/mbxha18/allcombinations/New_illumina/exeter/DAVE")
  list29903 <- read_table("29903.csv")
  all <- read_table("Copy of TWIST_synthetic_mixes_sample_sheet_150222.txt")
  #print(all)
  selectedRows <- all[grep(substr(d, 83, 85), all$plate_well), ]
  print(selectedRows)
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
  deltaC29control$lineage <- c("deltaAY2")
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
  files <- list.files(pattern = "\\.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv$")
  LRoutput <- read.table(files[1], header = TRUE)
  LRoutput <- separate(data = LRoutput, col = "Cons.Cov.Reads1.Reads2.Freq.P.value", into = c("Col", "Cov", "Reads1", "Reads2", "VarFreq", "P-value"), sep = "\\:")
  LRoutput$VarFreq <- as.numeric(sub("%", "", LRoutput$VarFreq, fixed = TRUE)) / 100
  #print(LRoutput)
  # Read the LR depths file
  depth <- list.files(pattern = "\\.ill.samtools.cov")
  LRdepths <- read.table(depth[1], col.names = c("NAME", "POS", "DEPTH"))
  
  # Define the list of controls
  control_list <- list(alphacontrol, betacontrol, deltaC23control, deltaC29control, omicroncontrol)
  allcontrols <- do.call(rbind, control_list)
  #print(allcontrols)
  # Compute the total expected frequency
  total_expected_frequency <- aggregate(cbind(expectedfrequency) ~ POS, data = allcontrols, sum)
  #print(total_expected_frequency)
  substrings <- c("alpha", "beta", "deltaC23", "deltaC29", "omicron")
  matches <- str_extract_all(List2, paste(substrings, collapse = "|"))[[1]]
  list3 <- unique(matches)
  # Filter the controls for the true SNPs
  truecontrols <- allcontrols[allcontrols$lineage %in% list3,]
  #print(truecontrols)  
  # Compute the TP and FP for the LR output
  LRvsTRUE <- merge(LRoutput, truecontrols, by.x = c("Position", "Var"), by.y = c("POS", "ALT"))
  TP_ILL <- as.data.frame(LRvsTRUE)
  assign(paste0("TP_ILL", substr(d, 11,100)), TP_ILL)
  
  
  FP_ILL <- anti_join(LRoutput, truecontrols, by = c("Position" = "POS"))
  assign(paste0("FP_ILL", substr(d, 11,100)), FP_ILL)
  # Compute the FN for the true SNPs
  FN_ILL <- anti_join(truecontrols, LRvsTRUE, by = c("POS" = "Position"))
  
  FN_ILL$DEPTH <- LRdepths$DEPTH[match(FN_ILL$POS, LRdepths$POS)]
  assign(paste0("FN_ILL", substr(d, 11,100)), FN_ILL)
  # Filter the FN based on coverage
  
  # Compute the TN for the LR output
  TNBEFORELR <- anti_join(list29903, truecontrols, by = c("POS" = "POS"))
  TN_ILL <- anti_join(TNBEFORELR, LRoutput, by = c("POS" = "Position"))
  assign(paste0("TN_ILL", substr(d, 11,100)), TN_ILL)
  # Add the observed frequency, depth, lineage, and sample columns to the total expected frequency data frame
  total_expected_frequency$observedfrequency <- LRoutput$VarFreq[match(total_expected_frequency$POS, LRoutput$Position)]
  total_expected_frequency$depth <- LRdepths$DEPTH[match(total_expected_frequency$POS, LRdepths$POS)]
  total_expected_frequency$lineage <- List2
  total_expected_frequency[is.na(total_expected_frequency)] <- 0
  total_expected_frequency$sample <- d
  
  # Append the total expected frequency data frame to the list
  total_ILL <- total_expected_frequency
  assign(paste0("total_ILL", substr(d, 11,100)), total_ILL)
}








all3ONT <- do.call(rbind, lapply( ls(patt="total_ONT"), get) )
all3ILL <- do.call(rbind, lapply( ls(patt="total_ILL"), get) )
all3ONT$sequencing <- "ONT"
all3ILL$sequencing <- "ILL"
all3all <- rbind(all3ONT,all3ILL)
library(ggplot2)

all3ONT <- all3ONT %>% filter(all3ONT$expectedfrequency > 0)
all3ILL <- all3ILL %>% filter(all3ILL$expectedfrequency > 0)

dave <- ggplot(data=all3ONT, aes(all3ILL$observedfrequency, all3ONT$observedfrequency)) + theme_classic() + xlim(c(0,1)) + ylim(c(0,1)) + ylab("ONT Observed Frequency") + xlab("Illumina Observed Frequency") + geom_point(aes(x= all3ONT$observedfrequency, y=all3ILL$observedfrequency)) #add alpha=all3ONT$depth if required
corr <- cor.test(x=all3ONT$observedfrequency, y=all3ILL$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr






ggsave("outputplease.png", dave)
