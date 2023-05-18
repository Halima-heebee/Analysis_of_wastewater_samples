#Loose Lab Code (MG), 17.05.2023. Locates false negative positions called using either/both technologies. Graph 1 displays the count of each false negative position by sequencing technology. Graph 2 displays the count of each false negative for each sequencing technology on a position-by-position basis.

dirs <- list.dirs("/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002", recursive = FALSE) #Provide directory containing relevant Illumina and ONT varscans. Directory structure should be one folder for each well (i.e. A1-NIMAGEN-001, A2-NIMAGEN-002 etc).
dirs <- dirs[!grepl("RECYCLE", dirs)] #removes the recycle bin when not running on terminal
library(readr) #loads appropriate packages into environment
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

for (d in dirs) {
  setwd("/home/mbxha18/allcombinations/New_illumina/exeter/DAVE") #Path containing list of positions in genome, Control SNP Locations, and the relationship between plate well and synthetic SARS-CoV-2 mixture
  list29903 <- read_table("29903.csv") #List of positions in genome
  all <- read_table("Copy of TWIST_synthetic_mixes_sample_sheet_150222.txt") #Matches plate well to expected frequencies of SARS-CoV-2 variants. Note - for exeter NT001, change this file to underscore.txt
  selectedRows <- all[grep(substr(d, 82, 84), all$plate_well), ] #Note - change these numbers dependent on the naming of the path - important for Exeter NT001.
  alphaEF <- as.numeric(selectedRows$alpha_exp_frequency[1]) #Extracts expected frequency from plate well data
  betaEF <- selectedRows$beta_exp_frequency[1]
  deltaC23EF <- selectedRows$deltaC23_exp_frequency[1]
  deltaC29EF <- selectedRows$deltaC29_exp_frequency[1]
  omicronEF <- selectedRows$omicron_exp_frequency[1]
  alphacontrol <- read_table("AlphaC15.tsv") #Loads Positions of Control SNPs
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
  expfrequencyonly <- expfrequencyonly[, colSums(expfrequencyonly != 0) > 0] #selects only SNPs with an expected frequency
  list <- colnames(expfrequencyonly)
  list <- gsub('_exp_frequency','', list)
  List2 <- paste(unlist(list),collapse="") #lists expected lineages for the sample
  # Set the working directory
  setwd(d)
  # Read the LR output file
  files <- list.files(pattern = "\\.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv$") #read in varscan - here you select which combination of thresholds you wish to use
  #print(files)
  LRoutput <- read.table(files[1], header=TRUE)
  LRoutput <- separate(data = LRoutput, col = "Cons.Cov.Reads1.Reads2.Freq.P.value", into = c("Col", "Cov", "Reads1", "Reads2", "VarFreq", "P-value"), sep = "\\:")
  LRoutput$VarFreq <- as.numeric(sub("%", "", LRoutput$VarFreq, fixed = TRUE)) / 100 #Removes % sign from VarFreq column of mpileup
  #print(LRoutput)
  # Read the LR depths file
  depth <- list.files(pattern = "\\.ont.samtools.cov") #read in coverage file, generated using SamTools
  LRdepths <- read.table(depth[1], col.names = c("NAME", "POS", "DEPTH"))
  
  # Define the list of controls
  control_list <- list(alphacontrol, betacontrol, deltaC23control, deltaC29control, omicroncontrol)
  allcontrols <- do.call(rbind, control_list)
  
  # Compute the total expected frequency
  total_expected_frequency <- aggregate(cbind(expectedfrequency) ~ POS, data = allcontrols, sum)
  substrings <- c("alpha", "beta", "deltaC23", "deltaC29", "omicron")
  matches <- str_extract_all(List2, paste(substrings, collapse = "|"))[[1]]
  list3 <- unique(matches)
  # Filter the controls for the true SNPs
  truecontrols <- allcontrols[allcontrols$lineage %in% list3,]
  
  # Compute the TP and FP for the LR output
  LRvsTRUE <- merge(LRoutput, truecontrols, by.x = c("Position", "Var"), by.y = c("POS", "ALT")) #identifies true positives
  TP_ONT <- as.data.frame(LRvsTRUE)
  assign(paste0("TP_ONT", substr(d, 11,100)), TP_ONT) #gives each sample (i.e. each plate well) a true positive dataframe
 
  
  FP_ONT <- anti_join(LRoutput, truecontrols, by = c("Position" = "POS")) #identifies false positives
  assign(paste0("FP_ONT", substr(d, 11,100)), FP_ONT)
  # Compute the FN for the true SNPs
  FN_ONT <- anti_join(truecontrols, LRvsTRUE, by = c("POS" = "Position")) #gives each sample (i.e. each plate well) a false negative dataframe
  
  FN_ONT$DEPTH <- LRdepths$DEPTH[match(FN_ONT$POS, LRdepths$POS)]
  assign(paste0("FN_ONT", substr(d, 11,100)), FN_ONT)
  # Filter the FN based on coverage
  #FNbecausecoverageONT <- FN_ONT[[d]][FN_ONT[[d]]$DEPTH < 3,] #can use to identify which false negatives arise from low coverage
  #FNbecausejustnotcalledONT <- FN_ONT[[d]][FN_ONT[[d]]$DEPTH > 3,]
  
  # Compute the TN for the LR output
  TNBEFORELR <- anti_join(list29903, truecontrols, by = c("POS" = "POS"))
  TN_ONT <- anti_join(TNBEFORELR, LRoutput, by = c("POS" = "Position")) #identifies true negatives
  assign(paste0("TN_ONT", substr(d, 11,100)), TN_ONT) #gives each sample (i.e. each plate well) a true negative dataframe

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
  setwd("/home/mbxha18/allcombinations/New_illumina/exeter/DAVE") #Does the exact same process for Illumina
  list29903 <- read_table("29903.csv")
  all <- read_table("Copy of TWIST_synthetic_mixes_sample_sheet_150222.txt") #Remember to change to underscore.txt for Exeter NT001
  selectedRows <- all[grep(substr(d, 82, 84), all$plate_well), ] #Careful
  alphaEF <- as.numeric(selectedRows$alpha_exp_frequency[1])
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
  expfrequencyonly <- expfrequencyonly[, colSums(expfrequencyonly != 0) > 0]
  list <- colnames(expfrequencyonly)
  list <- gsub('_exp_frequency','', list)
  List2 <- paste(unlist(list),collapse="")
  # Set the working directory
  setwd(d)
  # Read the LR output file
  files <- list.files(pattern = "\\.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv$") #Choose wisely
  #print(files)
  LRoutput <- read.table(files[1], header=TRUE)
  LRoutput <- separate(data = LRoutput, col = "Cons.Cov.Reads1.Reads2.Freq.P.value", into = c("Col", "Cov", "Reads1", "Reads2", "VarFreq", "P-value"), sep = "\\:")
  LRoutput$VarFreq <- as.numeric(sub("%", "", LRoutput$VarFreq, fixed = TRUE)) / 100
  #print(LRoutput)
  # Read the LR depths file
  depth <- list.files(pattern = "\\.ill.samtools.cov") #Manual input required here - a SamTools coverage file
  LRdepths <- read.table(depth[1], col.names = c("NAME", "POS", "DEPTH"))
  
  # Define the list of controls
  control_list <- list(alphacontrol, betacontrol, deltaC23control, deltaC29control, omicroncontrol)
  allcontrols <- do.call(rbind, control_list)
  
  # Compute the total expected frequency
  total_expected_frequency <- aggregate(cbind(expectedfrequency) ~ POS, data = allcontrols, sum)
  substrings <- c("alpha", "beta", "deltaC23", "deltaC29", "omicron")
  matches <- str_extract_all(List2, paste(substrings, collapse = "|"))[[1]]
  list3 <- unique(matches)
  # Filter the controls for the true SNPs
  truecontrols <- allcontrols[allcontrols$lineage %in% list3,]
  
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
  #FNbecausecoverageONT <- FN_ONT[[d]][FN_ONT[[d]]$DEPTH < 3,]
  #FNbecausejustnotcalledONT <- FN_ONT[[d]][FN_ONT[[d]]$DEPTH > 3,]
  
  # Compute the TN for the LR output
  TNBEFORELR <- anti_join(list29903, truecontrols, by = c("POS" = "POS"))
  TN_ILL <- anti_join(TNBEFORELR, LRoutput, by = c("POS" = "Position"))
  assign(paste0("TP_ILL", substr(d, 11,100)), TP_ILL)
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

#False Negatives

library(plyr)
df_list <- mget(ls(pattern = "FN_ILL.*")) #merge all Illuminas False Negative dataframes together
FN_ILL <- plyr::rbind.fill(df_list) 
FN_ILL$sequencing <- "ILL"
df_list <- mget(ls(pattern = "FN_ONT*")) #merge all ONTs False Negative dataframes together
FN_ONT <- plyr::rbind.fill(df_list)
FN_ONT$sequencing <- "ONT"
FN_ONT$count <- FN_ONT %>% add_count(POS) #add count of False Negatives at each position for that sequencing technology
FN_ILL$count <- FN_ILL %>% add_count(POS)

FN_ONT_distinct <- FN_ONT %>% 
  distinct(POS, count$n, sequencing, .keep_all = TRUE)
FN_ILL_distinct <- FN_ILL %>% 
  distinct(POS, count$n, sequencing, .keep_all = TRUE)

FN_ILL_distinct <- FN_ILL_distinct[order(FN_ILL_distinct$POS), ] #orders dataframe by position
FN_ONT_distinct <- FN_ONT_distinct[order(FN_ONT_distinct$POS), ]

FN_all_all <- full_join(FN_ONT_distinct, FN_ILL_distinct) #joins Illumina and ONT FNs

FN_all_all$seq_type <- ifelse(FN_all_all$sequencing=="ONT", "ONT", "ILL")

# Summarize the data by position and sequencing type
FN_all_all_summary <- aggregate(count$n ~ POS + seq_type, data=FN_all_all, sum)

# Create a new data frame to plot
df_plot <- data.frame(
  POS = FN_all_all_summary$POS,
  ONT = ifelse(FN_all_all_summary$seq_type=="ONT", FN_all_all_summary$`count$n`, 0),
  ILL = ifelse(FN_all_all_summary$seq_type=="ILL", FN_all_all_summary$`count$n`, 0)
)

df_plot <- df_plot %>% arrange(df_plot$POS)
duplicates <- df_plot$POS[!(duplicated(df_plot$POS)|duplicated(df_plot$POS, fromLast=TRUE))]
m <- as.data.frame(matrix(0,ncol = 3, nrow = length(duplicates)))
m$POS <- duplicates
m$ONT <- 0
m$ILL <- 0

m <- m %>% select(-c(V1, V2, V3))


everything_is_duplicated <- rbind(df_plot, m)
everything_is_duplicated <- everything_is_duplicated %>% arrange(everything_is_duplicated$ONT)
everything_is_duplicated <- everything_is_duplicated %>% arrange(everything_is_duplicated$POS)

shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}
everything_is_duplicated$ONT <- shift(everything_is_duplicated$ONT, 1)


toDelete <- seq(1, nrow(everything_is_duplicated), 2)
everything_is_duplicated <- everything_is_duplicated[ toDelete ,] #deletes duplicated False Negative position rows

ploteverything <- ggplot(everything_is_duplicated) +
  geom_point(aes(y = ONT, x = ILL)) +
  xlab("Illumina False Negative Count") +
  ylab("ONT False Negative Count") +
  ggtitle("Scatter plot of Illumina versus ONT False Negative counts") +
  theme_classic() + 
  geom_abline(linetype = "dotted") +
  geom_text(aes(x = ILL, y = ONT, label = ifelse(ONT > 20 | ILL > 20, as.character(POS), "")),
            vjust = -0.01, size = 3, hjust = -0.2,
            nudge_y = ifelse(everything_is_duplicated$POS == 29742, -3, 0)) + xlim(0,80) #Change which positions overlap using this line of code. Labels all positions with ONT and ILL False Negative Counts > 20.

ggsave("FalseNegatives_ExeterNT002.png",ploteverything) #Plot 


everything_is_duplicated <- everything_is_duplicated %>%
  filter(ONT >= 8 | ILL >= 8)

# Create the plot
positionbyposition <- ggplot(everything_is_duplicated) +
  geom_point(aes(y=ONT, x=ILL)) +
  xlab("Illumina False Negative Count") +
  ylab("ONT False Negative Count") +
  ggtitle("Scatter plot of Illumina versus ONT False Negative counts") +facet_wrap(~ POS, nrow = 5)

ggsave("FalseNegatives_ExeterNT002_PositionbyPosition.png",positionbyposition) #Plots counts of False Negatives using ONT and ILL technologies by position
