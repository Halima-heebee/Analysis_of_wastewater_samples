setwd("D:/")
install.packages("hrbrthemes")
library(readr)
library(ggplot2)
library(hrbrthemes)
library(dplyr)
library(tidyr)
library(viridis)
library(readxl)

# Read the table from Excel file
ILLtable <- read_excel("NOTTS_NT001_ILL.xlsx")
ONTtable <- read_excel("NOTTS_NT001_ONT.xlsx")


# Create a violin plot and boxplot for precision
p_meds <- ONTtable %>%
  group_by(threshold) %>%
  summarise(med = median(precision))

precision_plot <- ggplot(ONTtable, aes(x = threshold, y = precision, fill = threshold)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  ylim(0, 1) +
  geom_text(
    data = p_meds,
    aes(x = threshold, y = med, label = med),
    size = 2,
    vjust = -10
  )

# Save precision plot as PNG
ggsave("precision_plot.png", plot = precision_plot, width = 6, height = 4)
precision_plot
# Create a violin plot and boxplot for sensitivity
sensitivity_plot <- ggplot(ONTtable, aes(x = threshold, y = sensitivity, fill = threshold)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1) +
  theme_classic() +
  ylim(0, 1) +
  geom_text(
    data = p_meds,
    aes(x = threshold, y = med, label = med),
    size = 2,
    hjust = -0.5,
    vjust = -4
  )

# Save sensitivity plot as PNG
ggsave("sensitivity_plot.png", plot = sensitivity_plot, width = 6, height = 4)

# Filter the table for specific threshold values

ONTtable <- ONTtable %>% filter(threshold == "VF 0.01, MINCOV 1, MIN ALT READS 1, NO STRAND BIAS")
ILLtable <- ILLtable %>% filter(threshold == "VF 0.01, MINCOV 1, MIN ALT READS 1, NO STRAND BIAS")
bothtable <- rbind(ONTtable, ILLtable)

# Perform Wilcoxon signed-rank test for sensitivity
wilcox_result <- wilcox.test(ONTtable$sensitivity, ILLtable$sensitivity, paired = TRUE, exact = TRUE)

# Print the test result
print(wilcox_result)
