# Load necessary libraries
library(ggpubr)
library(dplyr)
library(tidyr)
# multi panel plot with patchwork
library(patchwork)

pyscomotif_fp1_data <- read.csv("result/SCOPe/250124_scope40_pyscomotif_benchmark_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)
folddisco_fp1_data <- read.csv("result/SCOPe/250215_scope_benchmark_with_lp_and_match.tsv.parsed.tsv", sep = "\t", header = TRUE)
# Join the two dataframes
pyscomotif_fp1_data$method <- "pyScoMotif"
folddisco_fp1_data$method <- "Folddisco"
data <- rbind(pyscomotif_fp1_data, folddisco_fp1_data)

tool_colors <- c('Folddisco' = "#FFA200", 'pyScoMotif' = "#D0D0D0")

recall_fp1_plot <- ggviolin(
    data[data$ratio %in% c(0.2,0.4,0.6,0.8,1.0), ], x = "ratio", y = "recall", fill = "method", 
    facet.by = "method",
    add = "mean_se",
    ylab = "Sensitivity at 1th FP",
    xlab = "Ratio",
    palette = tool_colors,
    ggtheme = theme_pubr(),
    legend = "none",
    ylim = c(0, 1)) + 
    stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), vjust = -0.25, hjust=-0.7, color = "black") +
    theme(strip.background = element_blank(), strip.text.x = element_blank())

recall_fp1_plot


recall_plot <- ggviolin(data[data$ratio %in% c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0), ], x = "ratio", y = "recall", fill = "darkorange",
                        add = "mean_se",
                            ylab = "Sensitivity at 5th FP",
                            xlab = "Ratio",
                            ggtheme = theme_pubr(),
                            ylim = c(0, 1)) + 
                        stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), vjust = -0.25, hjust=-0.7, color = "black")


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
ggsave("250118_recall_plot_fp5.pdf", plot = recall_plot, width = 8, height = 6, units = "in")

ggsave("250118_f1.png", plot = , width = 8, height = 6, units = "in")

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
