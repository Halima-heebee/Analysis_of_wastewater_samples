# Step 1: Read data from text file
#Concatenates all big tables awk 'FNR==1 && NR!=1{next;}{print $0 "\t" FILENAME}' bigtable*.tsv > bigoutput.tsv
# Step 1: Read data from text file
setwd("D:/")
data <- read.table("bigoutput.tsv", header = TRUE)

# Step 2: Extract relevant columns
data <- data[, c("sequencing", "Jaccard_similarity", "thresholds")]

# Step 3: Group data by sequencing platform and threshold
library(dplyr)
grouped_data <- data %>%
  group_by(sequencing, thresholds)

# Step 4: Calculate median Jaccard_similarity for each group
median_similarities <- aggregate(Jaccard_similarity ~ sequencing + thresholds, grouped_data, median)
# Step 5: Plot median Jaccard_similarity values
library(ggplot2)
ggplot(median_similarities, aes(x = thresholds, y = Jaccard_similarity, color = sequencing)) +
  geom_line() +
  labs(x = "Threshold", y = "Median Jaccard_similarity", color = "Sequencing platform") + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

