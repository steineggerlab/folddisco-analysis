# Load necessary libraries
library(ggpubr)
library(dplyr)
library(tidyr)
# multi panel plot with patchwork
library(patchwork)
# library(extrafont)
# Import and register system fonts
# font_import()
# loadfonts(device = "pdf")  # Load fonts for PDF device

benchmark_data <- read.csv("result/time/serine_query_accuracy.tsv", header = TRUE, sep = "\t")
benchmark_data$runtime <- c(0.529, 5.02, NA)

# Pivot data for easier plotting, excluding accuracy
benchmark_data_long <- benchmark_data %>%
  select(-accuracy) %>%
  pivot_longer(cols = precision:f1_score, names_to = "metric", values_to = "value")

# Specify order for the methods on the x-axis
benchmark_data_long$type <- factor(benchmark_data_long$type, levels = c("Folddisco", "pyScoMotif", "PDB"))

# Create the ggplot with a secondary y-axis for runtime
serine_query_accuracy <- ggbarplot(
    benchmark_data_long, x = "metric", y = "value",
    fill = "type", palette = c("#FFA200", "#D3D3D3", "#BEBEBE"), 
    label = TRUE, lab.pos = "in", 
    position = position_dodge(0.8),
    color = NA,  # Removes the bar borders
    ylab = "Metric Value", xlab = "",
    legend = "none" # No legend
    # Re-label the x-axis: Folddisco, pyScoMotif, PDB
) + theme(
    text = element_text(size = 18),
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
)
 
serine_query_runtime <- ggbarplot(
    benchmark_data[!is.na(benchmark_data$runtime),], x = "type", y = "runtime",
    fill = "type", palette = c("#FFA200", "#D3D3D3", "#BEBEBE"), 
    label = TRUE, lab.pos = "in", 
    color = NA,  # Removes bar borders
    ylab = "Runtime (s)", xlab = "",
    legend = "none" # No legend
) + theme(
    text = element_text(size = 18), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
)

merged <- serine_query_accuracy + serine_query_runtime + plot_layout(ncol = 2, widths = c(2, 1))
merged

# save the plot as pdf
ggsave("result/img/serine_query_accuracy.pdf", merged, bg = "transparent",  width = 12, height = 4)
