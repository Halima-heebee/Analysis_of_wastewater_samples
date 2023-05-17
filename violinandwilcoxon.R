setwd("D:/")
install.packages("hrbrthemes")
library(readr)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(readxl)
table <- read_excel("NOTTS_NT001_ONT.xlsx")

library(plyr)
p_meds <- ddply(table, .(threshold), summarise, med = median(precision))
ggplot(table, aes(x=threshold, y=precision, fill=threshold)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() + ylim(0,1) + geom_text(data = p_meds, aes(x = threshold, y = med, label = med), size = 2, vjust = -10) 
dave <- ggplot(table, aes(x=threshold, y=sensitivity, fill=threshold)) + 
  geom_violin(trim=FALSE) + geom_boxplot(width=0.1) + theme_classic() + ylim(0,1) + geom_text(data = p_meds, aes(x = threshold, y = med, label = med), size = 2, hjust = -0.5, vjust = -4) 

library(tidyverse)
library(ggpubr)
library(rstatix)
ggplot_build(dave)$data


#friedman.test(sensitivity ~ sample | Run, data=as.matrix(table))

ONTtable <- table %>% filter(threshold == "VF 0.01, MINCOV 1, MIN ALT READS 1, NO STRAND BIAS")
ILLtable <- table %>% filter(threshold == "VF 0.03, MINCOV 8, MIN ALT READS 6, STRAND BIAS")
bothtable <- rbind(ONTtable, ILLtable)
wilcox.test(ONTtable$sensitivity, ILLtable$sensitivity, paired=TRUE, exact=TRUE) 

