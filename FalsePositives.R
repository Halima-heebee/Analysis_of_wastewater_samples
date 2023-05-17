dirs <- list.dirs("/home/mbxha18/allcombinations/New_illumina/exeter/exeter_ill_vs_exeter_ont_NT002", recursive = FALSE)
dirs <- dirs[!grepl("RECYCLE", dirs)]

library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(ggrepel)
for (d in dirs) {
  setwd("/home/mbxha18/allcombinations/New_illumina/exeter/DAVE")
  list29903 <- read_table("29903.csv")
  all <- read_table("Copy of TWIST_synthetic_mixes_sample_sheet_150222.txt")
  selectedRows <- all[grep(substr(d, 82, 84), all$plate_well), ]
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
  files <- list.files(pattern = "\\.ont.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv$")
  #print(files)
  LRoutput <- read.table(files[1], header=TRUE)
  LRoutput <- separate(data = LRoutput, col = "Cons.Cov.Reads1.Reads2.Freq.P.value", into = c("Col", "Cov", "Reads1", "Reads2", "VarFreq", "P-value"), sep = "\\:")
  LRoutput$VarFreq <- as.numeric(sub("%", "", LRoutput$VarFreq, fixed = TRUE)) / 100
  #print(LRoutput)
  # Read the LR depths file
  depth <- list.files(pattern = "\\.ont.samtools.cov")
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
  LRvsTRUE <- merge(LRoutput, truecontrols, by.x = c("Position", "VarFreq"), by.y = c("POS", "ALT"))
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
  selectedRows <- all[grep(substr(d, 82, 84), all$plate_well), ]
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
  files <- list.files(pattern = "\\.ill.mincov1.mpileup2snp.nostrandbiasfilter.varscan.tsv$")
  #print(files)
  LRoutput <- read.table(files[1], header=TRUE)
  LRoutput <- separate(data = LRoutput, col = "Cons.Cov.Reads1.Reads2.Freq.P.value", into = c("Col", "Cov", "Reads1", "Reads2", "VarFreq", "P-value"), sep = "\\:")
  LRoutput$VarFreq <- as.numeric(sub("%", "", LRoutput$VarFreq, fixed = TRUE)) / 100
  #print(LRoutput)
  # Read the LR depths file
  depth <- list.files(pattern = "\\.ill.samtools.cov")
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
  LRvsTRUE <- merge(LRoutput, truecontrols, by.x = c("Position", "VarFreq"), by.y = c("POS", "ALT"))
  TP_ILL <- as.data.frame(LRvsTRUE)
  assign(paste0("TP_ILL", substr(d, 11,100)), TP_ILL)
  
  
  FP_ILL <- anti_join(LRoutput, truecontrols, by = c("Position" = "POS"))
  assign(paste0("FP_ILL", substr(d, 11,100)), FP_ILL)
  # Compute the FN for the true SNPs
  FN_ILL <- anti_join(truecontrols, LRvsTRUE, by = c("POS" = "Position"))
  
  FN_ILL$DEPTH <- LRdepths$DEPTH[match(FN_ILL$POS, LRdepths$POS)]
  assign(paste0("FN_ILL", substr(d, 11,100)), FN_ILL)
  # Filter the FN based on coverage
  #FNbecausecoverageILL <- FN_ILL[[d]][FN_ILL[[d]]$DEPTH < 3,]
  #FNbecausejustnotcalledILL <- FN_ILL[[d]][FN_ILL[[d]]$DEPTH > 3,]
  
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









#Observed versus Observed
StatClipEllipse <- ggproto("StatClipEllipse", Stat,
                           required_aes = c("x", "y"),
                           compute_group = function(data, scales, type = "t", level = 0.95,
                                                    segments = 51, na.rm = FALSE) {
                             xx <- ggplot2:::calculate_ellipse(data = data, vars = c("x", "y"), type = type,
                                                               level = level, segments = segments)
                             xx %>% mutate(x=pmax(x, 0)) %>%
                               mutate(x=pmin(x, 1))
                           }
)
stat_clip_ellipse <- function(mapping = NULL, data = NULL,
                              geom = "path", position = "identity",
                              ...,
                              type = "t",
                              level = 0.90,
                              segments = 51,
                              na.rm = FALSE,
                              show.legend = NA,
                              inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatClipEllipse,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      type = type,
      level = level,
      segments = segments,
      na.rm = na.rm,
      ...
    )
  )
}
ggplot(data=all3ONT, aes(all3ONT$observedfrequency, all3ILL$observedfrequency)) + theme_classic() + stat_clip_ellipse(fill="red", geom="polygon",level=0.95,alpha=0.6) + xlim(c(0,1)) + ylim(c(0,1)) + xlab("ONT Observed Frequency") + ylab("Illumina Observed Frequency") + geom_point(aes(x= all3ONT$observedfrequency, y=all3ILL$observedfrequency))
corr <- cor.test(x=all3ONT$observedfrequency, y=all3ILL$observedfrequency, method = 'pearson') #calculates Spearman's Rank Correlation Coefficient
corr

#False Positives

library(plyr)
df_list <- mget(ls(pattern = "FP_ILL.*"))
FP_ILL <- plyr::rbind.fill(df_list)
FP_ILL$sequencing <- "ILL"
df_list <- mget(ls(pattern = "FP_ONT*"))
FP_ONT <- plyr::rbind.fill(df_list)
FP_ONT$sequencing <- "ONT"
FP_ONT$count <- FP_ONT %>% add_count(Position)
FP_ILL$count <- FP_ILL %>% add_count(Position)

FP_ONT_distinct <- FP_ONT %>% 
  distinct(Position, count$n, sequencing, .keep_all = TRUE)
FP_ILL_distinct <- FP_ILL %>% 
  distinct(Position, count$n, sequencing, .keep_all = TRUE)

FP_ILL_distinct <- FP_ILL_distinct[order(FP_ILL_distinct$Position), ]
FP_ONT_distinct <- FP_ONT_distinct[order(FP_ONT_distinct$Position), ]

FP_all_all <- full_join(FP_ONT_distinct, FP_ILL_distinct)

FP_all_all$seq_type <- ifelse(FP_all_all$sequencing=="ONT", "ONT", "ILL")

# Summarize the data by position and sequencing type
FP_all_all_summary <- aggregate(count$n ~ Position + seq_type, data=FP_all_all, sum)

# Create a new data frame to plot
df_plot <- data.frame(
  Position = FP_all_all_summary$Position,
  ONT = ifelse(FP_all_all_summary$seq_type=="ONT", FP_all_all_summary$`count$n`, 0),
  ILL = ifelse(FP_all_all_summary$seq_type=="ILL", FP_all_all_summary$`count$n`, 0)
)

df_plot <- df_plot %>% arrange(df_plot$Position)
duplicates <- df_plot$Position[!(duplicated(df_plot$Position)|duplicated(df_plot$Position, fromLast=TRUE))]
m <- as.data.frame(matrix(0,ncol = 3, nrow = length(duplicates)))
m$Position <- duplicates
m$ONT <- 0
m$ILL <- 0

m <- m %>% select(-c(V1, V2, V3))


everything_is_duplicated <- rbind(df_plot, m)
everything_is_duplicated <- everything_is_duplicated %>% arrange(everything_is_duplicated$ONT)
everything_is_duplicated <- everything_is_duplicated %>% arrange(everything_is_duplicated$Position)

shift <- function(x, n){
  c(x[-(seq(n))], rep(NA, n))
}
everything_is_duplicated$ONT <- shift(everything_is_duplicated$ONT, 1)


toDelete <- seq(1, nrow(everything_is_duplicated), 2)
everything_is_duplicated <- everything_is_duplicated[ toDelete ,]
highcount <- everything_is_duplicated[everything_is_duplicated$ONT > 3 |  everything_is_duplicated$ILL > 3, ]    

# Create the plot
dave <- ggplot(highcount) +
  geom_point(aes(y=ONT, x=ILL)) +
  xlab("Illumina False Positive Count") +
  ylab("ONT False Positive Count") +
  ggtitle("Scatter Plot of ILL vs ONT counts") + 
  geom_abline(linetype = "dotted") + theme_classic()

ggsave("FalsePositives_ExeterNT002.png",dave)


dave1 <- ggplot(highcount) +
  geom_point(aes(y=ONT, x=ILL)) +
  xlab("Illumina False Positive Count") +
  ylab("ONT False Positive Count") +
  ggtitle("Scatter Plot of ILL vs ONT counts") +
  geom_abline(linetype = "dotted") +
  theme_classic() +
  ylim(0,100) + xlim(0,50) +
  geom_text(aes(x=ILL, y=ONT, label=ifelse(ONT>15 & ILL>15, as.character(Position), "")), size=3, hjust = 1.2)

ggsave("FalsePositives_ExeterNT002_labelled_positions.png",dave1)

dave1.2 <- ggplot(highcount) +
  geom_point(aes(y=ONT, x=ILL)) +
  xlab("Illumina False Positive Count") +
  ylab("ONT False Positive Count") +
  ggtitle("Scatter Plot of ILL vs ONT counts") +
  geom_abline(linetype = "dotted") +
  theme_classic() +
  ylim(5,15) + xlim(5,15) +
  geom_text(aes(x=ILL, y=ONT, label=ifelse(ONT>7 & ILL>7, as.character(Position), "")), size=3, hjust = 1.2, position = position_jitter(width = 0.3, height = 0.3))
ggsave("FalsePositives_ExeterNT002_labelled_positions_jitter.png", dave1.2)


dave2 <- ggplot(highcount) +
  geom_point(aes(y=ONT, x=ILL)) +
  xlab("Illumina False Positive Count") +
  ylab("ONT False Positive Count") +
  ggtitle("Scatter Plot of ILL vs ONT counts") +
  geom_abline(linetype = "dotted") +
  theme_classic() +
  ylim(0,100) + xlim(0,50) +
  geom_text(data = subset(highcount, ILL >8 & ONT <= 5),
            aes(x = ILL, y = ONT, label = ifelse(Position == 18304, "18304", as.character(Position)),
                hjust = ifelse(Position == 18304, -0.7, -0.2)),
            size = 2.5, vjust = -0.5)

ggsave("FalsePositives_ExeterNT002_Illumina_labelled.png",dave2)

dave3 <- ggplot(highcount) +
  geom_point(aes(y=ONT, x=ILL)) +
  xlab("Illumina False Positive Count") +
  ylab("ONT False Positive Count") +
  ggtitle("Scatter Plot of ILL vs ONT counts") +
  geom_abline(linetype = "dotted") +
  theme_classic() +
  ylim(0,100) + xlim(-5,50) +
  geom_text(data = subset(highcount, ONT > 50 & ILL <= 0),
            aes(x = ILL, y = ONT, label = as.character(Position)),
            size = 1.5, hjust = 1.5, vjust = 0.2, check_overlap = TRUE,
            position = position_nudge(x = ifelse(duplicated(highcount[, c("ILL", "ONT")]), 2, 0))) +
  scale_x_continuous(breaks = seq(0, 50, by = 10))

ggsave("FalsePositives_ExeterNT002_ONT_Labelled.png", dave3)

highcount <- highcount %>%
  filter(ONT >= 30 | ILL >= 8)

dave4 <- ggplot(highcount) +
  geom_point(aes(y=ONT, x=ILL)) +
  xlab("Illumina False Positive Count") +
  ylab("ONT False Positive Count") +
  ggtitle("Scatter Plot of ILL vs ONT counts") + facet_wrap(~ Position, nrow = 5) +
  theme(axis.text.x = element_text(size=6, angle = 90, vjust = 0.5))
ggsave("FalsePositives_ExeterNT002_PositionbyPosition.png",dave4)


