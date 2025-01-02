# Load necessary libraries
library(ggpubr)
library(dplyr)
library(tidyr)
# multi panel plot with patchwork
library(patchwork)

data <- read.csv("241105_scope40_folddisco_benchmark_result_fp5.tsv.parsed.tsv", sep = "\t", header = TRUE)
data <- read.csv("241105_scope40_folddisco_benchmark_result_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)
# data <- read.csv("scope40_folddisco_benchmark_result.tsv.parsed.tsv", sep = "\t", header = TRUE)


# Plot precision 
# Boxplot: X-axis: ratio, y-axis: precision
# precision_plot <- ggboxplot(data, x = "ratio", y = "precision", 
#                             add = "jitter", 
#                             ylab = "Precision",
#                             xlab = "Ratio",
#                             ggtheme = theme_pubr())
recall_plot <- ggviolin(data[data$ratio %in% c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), ], x = "ratio", y = "recall", fill = "darkorange",
                        add = "mean_se",
                            ylab = "Sensitivity at 5th FP",
                            xlab = "Ratio",
                            ggtheme = theme_pubr(),
                            ylim = c(0, 1)) + 
                        stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), vjust = -0.25, hjust=-0.7, color = "black")
recall_plot

data$tpratio <- data$TP / data$answer_length


precision_plot <- ggviolin(data, x = "ratio", y = "precision", fill = "darkorange",
                        add = "mean_se",    
                            ylab = "Precision at 5th FP",
                            xlab = "Ratio",
                            ggtheme = theme_pubr(),
                            ylim = c(0, 1))
precision_plot

precision_plot + recall_plot

ggsave("precision_recall_plot.pdf", plot = precision_plot + recall_plot, width = 10, height = 5, units = "in")
ggsave("scope_benchmark/recall_plot_fp5.png", plot = recall_plot, width = 8, height = 6, units = "in")
# # Print where recall is over 0.3
# data %>% filter(recall > 0.6) %>% print


# recall_plot <- ggboxplot(data,
#     x = "ratio", y = "recall",
#     add = "jitter",
#     ylab = "Recall at 5th FP",
#     xlab = "Ratio",
#     ggtheme = theme_pubr(),
#     add.params = list(alpha = 0.15, color = "darkgray"),
#     ylim = c(0, 1)
# )
# recall_plot

# precision_plot <- ggboxplot(data,
#     x = "ratio", y = "precision",
#     add = "jitter",
#     ylab = "Precision at 5th FP",
#     xlab = "Ratio",
#     ggtheme = theme_pubr(),
#     add.params = list(alpha = 0.15, color = "darkgray"),
#     ylim = c(0, 1)
# )
# precision_plot

library(ggridges)

data$ratio_factor <- as.factor(rev(data$ratio))
data$recall_factor <- as.factor(data$recall)
recall_plot <- ggplot(data[data$answer_length > 3,], aes(x = recall, y = ratio_factor)) +
    geom_density_ridges(fill = "steelblue", scale = 1.5, quantile_lines = TRUE) + 
    # scale_x_continuous(limits = c(0, 1)) +
    theme_pubr()
recall_plot
