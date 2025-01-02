# Load necessary libraries
library(ggpubr)
library(dplyr)
library(tidyr)
# multi panel plot with patchwork
library(patchwork)

zinc_data <- read.csv("result/time/zinc_query_accuracy.tsv", header = TRUE, sep = "\t")
zinc_data$type <- factor(zinc_data$type, levels = c("Folddisco", "pyScoMotif", "PDB"))


# Filter data for query_len == 3
zinc_data_3 <- zinc_data %>% filter(query_len == 3)

# Filter data for query_len == 4
zinc_data_4 <- zinc_data %>% filter(query_len == 4)

# Pivot data for metrics plot (precision, recall, f1_score) for query_len == 3
zinc_data_long_3 <- zinc_data_3 %>%
  select(-accuracy) %>%
  pivot_longer(cols = precision:f1_score, names_to = "metric", values_to = "value")

# Pivot data for metrics plot (precision, recall, f1_score) for query_len == 4
zinc_data_long_4 <- zinc_data_4 %>%
  select(-accuracy) %>%
  pivot_longer(cols = precision:f1_score, names_to = "metric", values_to = "value")

# zinc_data_long_3$type <- factor(zinc_data_long_3$type, levels = c("Folddisco", "pyScoMotif", "PDB"))
# zinc_data_long_4$type <- factor(zinc_data_long_4$type, levels = c("Folddisco", "pyScoMotif", "PDB"))
# Create the metrics plot for query_len == 3 (Precision, Recall, F1-score)
zinc_plot_3_metrics <- ggbarplot(
  zinc_data_long_3, x = "metric", y = "value", 
  fill = "type", palette = c("#FFA200", "#D3D3D3", "#BEBEBE"), 
  label = TRUE, lab.pos = "in", 
  position = position_dodge(0.8),
  color = NA,  # Removes the bar borders
  ylab = "Metric Value", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = 18),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
)

# Create the metrics plot for query_len == 4 (Precision, Recall, F1-score)
zinc_plot_4_metrics <- ggbarplot(
  zinc_data_long_4, x = "metric", y = "value", 
  fill = "type", palette = c("#FFA200", "#D3D3D3", "#BEBEBE"), 
  label = TRUE, lab.pos = "in", 
  position = position_dodge(0.8),
  color = NA,  # Removes the bar borders
  ylab = "Metric Value", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = 18),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
)

# Create runtime plot for query_len == 3
zinc_plot_3_runtime <- ggbarplot(
  zinc_data_3[!is.na(zinc_data_3$runtime),], x = "type", y = "runtime",
  fill = "type", palette = c("#FFA200", "#D3D3D3", "#BEBEBE"), 
  label = TRUE, lab.pos = "in", 
  color = NA,  # Removes bar borders
  ylab = "Runtime (s)", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = 18),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
)

# Create runtime plot for query_len == 4
zinc_plot_4_runtime <- ggbarplot(
  zinc_data_4[!is.na(zinc_data_4$runtime),], x = "type", y = "runtime",
  fill = "type", palette = c("#FFA200", "#D3D3D3", "#BEBEBE"), 
  label = TRUE, lab.pos = "in", 
  color = NA,  # Removes bar borders
  ylab = "Runtime (s)", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = 18),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
)

# Combine plots for query_len == 3
merged_3 <- zinc_plot_3_metrics + zinc_plot_3_runtime + plot_layout(ncol = 2, widths = c(2, 1))

# Combine plots for query_len == 4
merged_4 <- zinc_plot_4_metrics + zinc_plot_4_runtime + plot_layout(ncol = 2, widths = c(2, 1))

# Display the plots
merged_3
merged_4

# Save the plots as PDF with transparent background
ggsave("result/img/zinc_finger_query_len_3.pdf", merged_3, bg = "transparent", width = 12, height = 4)
ggsave("result/img/zinc_finger_query_len_4.pdf", merged_4, bg = "transparent", width = 12, height = 4)
