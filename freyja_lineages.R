# load the ggplot2 library
library(ggplot2)
setwd("D:/")
# create a data frame with the summarized rows and lineage abundance
library("readr")
data <- read_csv("strainsinput.csv")
setwd("D:/")
library(readxl)
all <- read_excel("freyjatest.xlsx")
all$b.1.1.7_exp_frequency <- as.numeric(all$b.1.1.7_exp_frequency)
all$ba.1_exp_frequency <- as.numeric(all$ba.1_exp_frequency)

# reshape the data frame to a long format
data_merged <- merge(data, all, by = "plate_well")
data_merged <- data_merged[seq(1,nrow(data_merged), by=2),]
data_merged$stack <- tolower(data_merged$stack)
case4 <- data_merged$stack == "ba.1"
data_merged$freq <- ifelse(data_merged$stack == "b.1.1.17", data_merged$b.1.1.7_exp_frequency, data_merged$value)
library(dplyr)
new_df <- data_merged %>% mutate(frequency = case_when(stack == "b.1.1.7" ~ b.1.1.7_exp_frequency, stack == "b.1.351" ~ b.1.351_exp_frequency, stack == "b.1.617.2" ~ b.1.617.2_exp_frequency, stack == "ay.2" ~ ay.2_exp_frequency, case4 ~ ba.1_exp_frequency, TRUE ~ 0))
exp <- distinct(new_df, stack, plate_well, .keep_all = TRUE)
exp$group <- "EXP"
exp <- exp %>% mutate(value=frequency)
exp <- exp %>% select(-frequency)
final <- rbind(data_merged, exp)
final <- na.omit(final)

library(ggplot2)
library(stringr)
final$sample_label <- str_replace(final$SAMPLE, ".*([A-Za-z0-9]{3})\\.variants\\.tsv", "\\1")
library(dplyr)
final <- final %>% 
  group_by(group, sample_label) %>% 
  mutate(value = ifelse(is.na(value), 0, value), # replace NA values with 0
         stack = factor(stack, levels=unique(stack)), # ensure the stack factor levels are ordered consistently
         stack = reorder(stack, -value)) # reorder the levels within each group based on the value of value





#lighter means child lineage

final2 <- final
final2$color <- "DAVE"
final2$color <- ifelse(grepl("q.",final2$stack), "#F60C0C", final2$color)
final2$color <- ifelse(grepl("q.",final2$stack), "#F60C0C", final2$color)
final2$color <- ifelse(grepl("b.1.1.7", final2$stack), "#F8766D", final2$color)
final2$color <- ifelse(grepl("b.1.343", final2$stack), "#CD9500", final2$color)
final2$color <- ifelse(grepl("b.1.351", final2$stack), "#CD9600", final2$color)
final2$color <- ifelse(grepl("b.1.617.2", final2$stack), "#7CAE00", final2$color)
final2$color <- ifelse(grepl("ay.",final2$stack), "#1AA727", final2$color)
final2$color <- ifelse(grepl("ay.2",final2$stack), "#0E5615", final2$color)
final2$color <- ifelse(grepl("b.1.525",final2$stack), "#00BF7D", final2$color)
final2$color <- ifelse(grepl("p.1", final2$stack), "#00BFC4", final2$color)
final2$color <- ifelse(grepl("c.37", final2$stack), "#00A9FF", final2$color)
final2$color <- ifelse(grepl("c.38", final2$stack), "#00A9FF", final2$color)
final2$color <- ifelse(grepl("ba.2.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bc.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bd.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bg.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bh.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bj.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bl.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bm.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bn.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bp.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("br.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("bs.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("by.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("ca.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("cb.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("ch.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("cj.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("cm.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("cv.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("dd.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("ds.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("dv.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("ej.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("ep.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("fj.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("fk.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("fr.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("fs.",final2$stack), "#E76BF3", final2$color)
final2$color <- ifelse(grepl("ba.1", final2$stack), "#A58AFF", final2$color)
final2$color <- ifelse(grepl("DAVE", final2$color), "#FF61CC", final2$color)


final2$color <- as.factor(final2$color)
#color_levels <- levels(final2$color)

# Specify the desired order of the factor levels
#new_color_levels <- c("#F8766D", "#F60C0C", "#CD9600", "#7CAE00", "#1AA727", "#0E5615", "#A58AFF", "#FF61CC")

# Reorder the factor levels so that they match the desired order
#final2$color <- factor(final2$color, levels = new_color_levels)

#final2 <- transform(final2, color=factor(color, levels=unique(color)))

p1 <- ggplot(final2, aes(x=group, y=value, fill=color)) + 
  geom_bar(stat="identity", position="stack") + 
  facet_grid(. ~ sample_label) + 
  theme_classic()  + labs(x = "Sequencing Technology", y = "Relative Abundance") + scale_fill_identity(name="color", guide = 'legend', labels=c("#F8766D" = "alpha", "#F60C0C" = "children of alpha", "#CD9600" = "beta", '#7CAE00'='deltaC23', '#1AA727' = "children of deltaC23", "#0E5615" = "deltaC29", "#A58AFF" = "omicron", "#FF61CC" = "others"))


p1

p1
setwd("D:/")
ggsave("delta_omicron_1_L-delta_omicron_7_L.png", plot = p1, width = 10, height = 5, dpi = 300)