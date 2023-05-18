#Loose Lab Code (MG), 17.05.2023. Graph 1 displays the effect of the minimum concentration of a particular lineage on the sensitivity to SNPs of that lineage for each sequencing technology. Graph 2 displays whether pure or mixed sample wells provide optimum sensitivity values for each technology. Graph 3 displays the effect of both lineage and mixed nature of the sample.
dirs <- list.dirs("/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT001", recursive = FALSE) #Provide path to directories containing both ONT and Illumina varscans for each sample well
dirs <- dirs[!grepl("RECYCLE", dirs)] #removes the recycle bin when not running on terminal
library(readr) #loads appropriate packages into environment
library(dplyr)
library(stringr)
library(tidyr)
for (d in dirs) {
  setwd("/home/mbxha18/allcombinations/New_illumina/exeter/DAVE") #Path containing list of positions in genome, Control SNP Locations, and the relationship between plate well and synthetic SARS-CoV-2 mixture
  list29903 <- read_table("29903.csv") #List of positions in genome
  all <- read_table("underscore.txt") #Matches plate well to expected frequencies of SARS-CoV-2 variants. Note - for exeter NT002 and 3, change this file to Copy of TWIST_synthetic_mixes_sample_sheet_150222.txt.
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
  deltaC29control$lineage <- c("deltaC29")
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
  files <- list.files(pattern = "\\.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv$") #read in varscan - here you select which combination of thresholds you wish to use
  #print(files)
  LRoutput <- read.table(files[1], header=TRUE)
  LRoutput <- separate(data = LRoutput, col = "Cons.Cov.Reads1.Reads2.Freq.P.value", into = c("Col", "Cov", "Reads1", "Reads2", "VarFreq", "P-value"), sep = "\\:")
  LRoutput$VarFreq <- as.numeric(sub("%", "", LRoutput$VarFreq, fixed = TRUE)) / 100
  #print(LRoutput)
  # Read the LR depths file
  depth <- list.files(pattern = "\\.ill.samtools.cov") #read in coverage file, generated using SamTools
  LRdepths <- read.table(depth[1], col.names = c("NAME", "POS", "DEPTH"))
  
  
  # Define the list of controls
  control_list <- list(alphacontrol, betacontrol, deltaC23control, deltaC29control, omicroncontrol)
  allcontrols <- do.call(rbind, control_list)
  
  # Compute the total expected frequency
  total_expected_frequency <- aggregate(cbind(expectedfrequency) ~ POS, data = allcontrols, sum)
  substrings <- c("alpha", "beta", "deltaC23", "deltaC29", "omicron")
  matches <- str_extract_all(List2, paste(substrings, collapse = "|"))[[1]]
  list3 <- unique(matches)
  for (item in list3){
    truecontrols <- allcontrols[allcontrols$lineage == item, ]
    LRvsTRUE <- merge(LRoutput, truecontrols, by.x = c("Position", "Var"), by.y = c("POS", "ALT")) #identifies true positives
    TP_ONT <- as.data.frame(LRvsTRUE)
    assign(paste0("TP_ONT", substr(d, 11,100), "_", item), TP_ONT)
  # Compute the FN for the true SNPs
    FN_ONT <- anti_join(truecontrols, LRvsTRUE, by = c("POS" = "Position"))
    FN_ONT$DEPTH <- LRdepths$DEPTH[match(FN_ONT$POS, LRdepths$POS)]
    assign(paste0("FN_ONT", substr(d, 11,100), "_", item), FN_ONT)
  #addlineagename to the end of FP_ONT
     sensitivity <- nrow(TP_ONT)/(nrow(TP_ONT) + nrow(FN_ONT)) #calculates sensitivity for the sample
     assign(paste0("sensitivity", substr(d, 11,100), "_", item), sensitivity)
     sensitivity_df <- data.frame(matrix(ncol=3, nrow = 1))
     colnames(sensitivity_df) <- c("sensitivity", "min_conc", "lineage")
     sensitivity_df$sensitivity <- sensitivity
     list4 <- selectedRows[[item]]
     list5 <- item
     sensitivity_df$min_conc <- list4[1] #make the variable stored in item the column in selectedRows it searches for, not the actual string item
     sensitivity_df$lineage <- list5
     assign(paste0("sensitivityDF", substr(d, 11,100), "_", item), sensitivity_df)
     }
}

library(plyr)
all_objects <- ls()
target_objects <- grep("*sensitivityDF", all_objects, value=TRUE) #extracts all sensitivity dataframes. Can be modified for precision and Jaccard_similarity.

target_objects
df_list <- mget(target_objects)
final_df <- do.call(rbind,df_list)
library(ggplot2)

sensitivityvsminconc <- ggplot(data=final_df) + 
  geom_point(aes(x=min_conc, y=sensitivity, color=lineage)) +
  geom_smooth(aes(x=min_conc, y=sensitivity, color=lineage), method = "loess", se=FALSE) +
  scale_x_continuous(limits = c(0, 200)) +
  scale_color_manual(values = c("#F8766D", "#CD9600", "#7CAE00", "#0E5615", "#A58AFF")) +
  xlab("Minimum Concentration (Genome Copies per Microlitre)") + ylab("Sensitivity") +
  theme_classic() #only looks at samples with maximum 200 genome copies per microlitres. #Colors designed to match Freyja.

ggsave("Sensitivity_by_Lineage_ILL_ExeterNT001_Lineage.png",sensitivityvsminconc)


my_list <- Map(cbind, df_list, sample = names(df_list))
new_final_df <- do.call(rbind,my_list)
new_final_df$sample <- substr(new_final_df$sample, 85,87) #Change these numbers dependent on your data.
new_final_df$combo <- 0
new_final_df$combo<- ifelse(duplicated(new_final_df$sample) | duplicated(new_final_df$sample, fromLast=TRUE), "MIXED", "PURE") #if there are multiple lineages with expected frequencies > 0 for a sample, set the value of this column equal to MIXED, if not set the value equal to PURE
print(new_final_df)
sensiivityboth <- ggplot(data=new_final_df) + 
  geom_point(aes(x=min_conc, y=sensitivity, color=lineage, shape=combo)) +
  geom_smooth(aes(x=min_conc, y=sensitivity, color=lineage), method = "loess", se=FALSE) +
  scale_x_continuous(limits = c(0, 200)) +
  scale_color_manual(values = c("#F8766D", "#CD9600", "#7CAE00", "#0E5615", "#A58AFF")) +
  xlab("Minimum Concentration (Genome Copies per Microlitre)") + ylab("Sensitivity") +
  theme_classic() #only looks at samples with maximum 200 genome copies per microlitres

ggsave("Sensitivity_by_Lineage_ILL_ExeterNT001_Lineage_and_Combo.png",sensitivityboth)


sensitivitycombo <- ggplot(data=new_final_df) + 
  geom_point(aes(x=min_conc, y=sensitivity, color=combo, shape=combo)) +
  geom_smooth(aes(x=min_conc, y=sensitivity, color=combo), method = "loess", se=FALSE) +
  scale_x_continuous(limits = c(0, 200)) +
  scale_color_manual(values = c("#FF0000", "#0000FF")) +
  xlab("Minimum Concentration (Genome Copies per Microlitre)") + ylab("Sensitivity") +
  theme_classic() #only looks at samples with maximum 200 genome copies per microlitres

ggsave("Sensitivity_by_Lineage_ILL_ExeterNT001_Combo.png",sensitivitycombo)
