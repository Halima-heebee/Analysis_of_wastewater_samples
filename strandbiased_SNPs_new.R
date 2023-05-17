setwd("D:/")
library("readr")
library("ggplot2")

strand_biased_SNPS <- read_table("biasedSNPs_counts.tsv")
strand_biased_SNPS$Position <- as.numeric(strand_biased_SNPS$Position)
strand_biased_SNPS <- strand_biased_SNPS[complete.cases(strand_biased_SNPS), ]
#strand_biased_SNPS <- subset(strand_biased_SNPS, Count > 20)
library(ggplot2)
ggplot(strand_biased_SNPS, aes(x = Position, y = Count)) + 
  geom_bar(stat = "identity", width = 20) +
  scale_x_continuous(limits = c(1, 30000), breaks = seq(0, 30000, by = 2000)) +
  labs(x = "Position", y = "Count")  + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  geom_text(aes(label = ifelse(Count > 45, as.character(Position), "")), vjust = -0.5, size = 2.5, hjust = -0.05, 
            nudge_x = ifelse(strand_biased_SNPS$Position == 3370, 800, 
                             ifelse(strand_biased_SNPS$Position == 16428, -1100, 0)))

