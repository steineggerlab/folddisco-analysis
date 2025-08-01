# 00. Libraries
# Load necessary libraries
library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# multi panel plot with patchwork
library(patchwork)
folddisco_fp1_data <- read.csv("benchmark/result/SCOPe/250214_scope_benchmark_with_length_penalty_0.7_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)
folddisco_fp1_data <- read.csv("benchmark/result/SCOPe/250125_scope40_folddisco_benchmark_result_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)
font_size <- 7
ggviolin(
    folddisco_fp1_data[folddisco_fp1_data$ratio %in% c(0.2,0.4,0.6,0.8,1.0), ], x = "ratio", y = "recall",  
    add = "mean",
    ylab = "Sensitivity at 1th FP",
    xlab = "Ratio",
    # palette = tool_colors,
    ggtheme = theme_pubr(),
    legend = "none",
    size = 0.3,
    add.params = list(size = 0.2, shape = 20),
    font.label = list(size = 1.6),
    ylim = c(0, 1)) + 
    stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2), size = 0.5), vjust = -0.25, hjust=-0.7, color = "black") +
    stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.5, color = "black") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
      legend.position = "none",  # Removes the legend
      plot.title = element_blank(),  # Removes the title
      text = element_text(size = font_size), # Set all font size to 24
      panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
      plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
      legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
      legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
      strip.background = element_blank(), strip.text.x = element_blank() # No facet labels
    ) 
