# Load necessary libraries
library(ggpubr)
library(dplyr)
library(tidyr)
# multi panel plot with patchwork
library(patchwork)

data <- read.csv("temp/temp_result2.txt", sep="\t", header=TRUE)

# Plot precision 
# Boxplot: X-axis: ratio, y-axis: precision
# precision_plot <- ggboxplot(data, x = "ratio", y = "precision", 
#                             add = "jitter", 
#                             ylab = "Precision",
#                             xlab = "Ratio",
#                             ggtheme = theme_pubr())
recall_plot <- ggboxplot(data, x = "ratio", y = "recall", 
                            add = "jitter", 
                            ylab = "Recall at 5th FP",
                            xlab = "Ratio",
                            ggtheme = theme_pubr(),
                            add.params= list(alpha = 0.15, color = "darkgray"),
                            ylim = c(0, 1))
recall_plot

precision_plot <- ggboxplot(data, x = "ratio", y = "precision", 
                            add = "jitter", 
                            ylab = "Precision at 5th FP",
                            xlab = "Ratio",
                            ggtheme = theme_pubr(),
                            add.params= list(alpha = 0.15, color = "darkgray"),
                            ylim = c(0, 1))
precision_plot

precision_plot + recall_plot

ggsave("precision_recall_plot.pdf", plot = precision_plot + recall_plot, width = 10, height = 5, units = "in")

# # Print where recall is over 0.3
# data %>% filter(recall > 0.6) %>% print
