# 00. Libraries
# Load necessary libraries
library(ggpubr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# multi panel plot with patchwork
library(patchwork)

# 01. Setup
# 01-1. Load Data

# index benchmark for Figure 1b, 1c, and 1e
index_benchmark_data <- read.csv("result/indexing/index_folddisco.tsv", header = TRUE, sep = "\t")
additional_data <- read.csv("result/indexing/index_pyscomotif_pdb.tsv", header = TRUE, sep = "\t")
# Ensure mode is a factor
index_benchmark_data$mode <- factor(index_benchmark_data$index_mode, levels = c("normal", "big"))
index_benchmark_data$tool <- "Folddisco"
# Filtering & combinding index benchmark data for Figure 1b, 1c, and 1e
# Filter for datasets with 12 threads
filtered_data <- index_benchmark_data %>%
  filter(description %in% c('default', 'other model', 'database') & num_threads == 12) %>%
  arrange(num_structures)
# Combine the data
filtered_data <- bind_rows(filtered_data, additional_data)
# Ensure tool is a factor
filtered_data$tool <- factor(filtered_data$tool, levels = c('Folddisco', 'pyScoMotif', 'RCSB'))
# Create a mapping for the custom x-axis labels
filtered_data$data <- factor(filtered_data$data, 
                             levels = c('e_coli', 's_cerevisiae', 'h_sapiens', 'pdb', 'swissprot'),
                             labels = c('E. coli', 'Yeast', 'Human', 'PDB', 'Swissprot'))

# foldcomp data for Figure 1d
foldcomp_data <- read.csv("result/indexing/foldcomp_comparison.tsv", header = TRUE, sep = "\t")
# Calculate the total index size in GB
foldcomp_data <- foldcomp_data %>%
  group_by(data) %>%
  mutate(total_size_in_gb = sum(size_in_gb))
foldcomp_data.tool <- factor(foldcomp_data$tool, levels = c('Folddisco', 'pyScoMotif'))

# query benchmark for Figure 1e - 1l? 
query_benchmark_data <- read.csv("result/querying/query_folddisco.tsv", header = TRUE, sep = "\t")
# Pivot data for easier plotting, excluding accuracy
query_benchmark_data_long <- query_benchmark_data %>%
  select(-accuracy) %>%
  pivot_longer(cols = precision:f1_score, names_to = "metric", values_to = "value")

# Split data for each figure
serine_data <- query_benchmark_data %>% filter(rank == "S01")
serine_data$type <- factor(serine_data$type, levels = c("Folddisco", "pyScoMotif", "RCSB"))
zinc_data_3 <- query_benchmark_data %>% filter(query_len == 3 & rank == "C2H2")
zinc_data_3$type <- factor(zinc_data_3$type, levels = c("Folddisco", "pyScoMotif", "RCSB"))
zinc_data_4 <- query_benchmark_data %>% filter(query_len == 4 & rank == "C2H2")
zinc_data_4$type <- factor(zinc_data_4$type, levels = c("Folddisco", "pyScoMotif", "RCSB"))
serine_data_long <- query_benchmark_data_long %>% filter(rank == "S01")
serine_data_long$type <- factor(serine_data_long$type, levels = c("Folddisco", "pyScoMotif", "RCSB"))
zinc_data_3_long <- query_benchmark_data_long %>% filter(query_len == 3 & rank == "C2H2")
zinc_data_3_long$type <- factor(zinc_data_3_long$type, levels = c("Folddisco", "pyScoMotif", "RCSB"))
zinc_data_4_long <- query_benchmark_data_long %>% filter(query_len == 4 & rank == "C2H2")
zinc_data_4_long$type <- factor(zinc_data_4_long$type, levels = c("Folddisco", "pyScoMotif", "RCSB"))


# 01-2. Set fill colors for each tool
tool_colors <- c('Folddisco' = "#FFA200", 'pyScoMotif' = "#D3D3D3", 'RCSB' = "#AAAAAA")
font_size <- 7
label_size <- 1.6

# 02. Plotting
# Figure 1b
# Create the plot for number of structures
# ggtexttable
fig_1b_num_structures <- ggtexttable(
  filtered_data %>% distinct(data, num_structures),
  rows = NULL, cols = c("DB", "Structures"),
  theme = ttheme(
    base_style = "blank", base_size = font_size,
    padding = unit(c(0.8, 2), "mm")
)) %>% tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)

# Figure 1c
# Create the plot for Index Size
fig_1c_index_size <- ggline(
  filtered_data, x = "data", y = "total_index_size_in_gb", color = "tool",
  palette = tool_colors,  # Set color for normal and big mode
  label = TRUE, lab.pos = "out", lab.size = label_size,  # Move annotations to the top of bars
  label.nb.digits = 1,
) +
  labs(
    y = "Index Size (GB)", x = ""
  ) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = font_size), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  ) + coord_cartesian(clip = 'off')


# Figure 1d
# Create the plot as stacked bar plot
fig_1d_foldcomp <- ggplot(foldcomp_data, aes(x = tool, y = size_in_gb, fill = tool)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  scale_fill_manual(values = tool_colors) +
  labs(y = "Size in GB", x = "") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = font_size), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  ) # Change the x-axis labels
fig_1d_foldcomp <- fig_1d_foldcomp + scale_x_discrete(labels = c("Ours", "pyScoMotif"))

# Create the plot for Runtime
fig_1e_runtime <- ggline(
  filtered_data %>% mutate(runtime_in_sec = round(runtime_in_sec, 0)),
  x = "data", y = "runtime_in_sec", color = "tool", 
  palette = tool_colors,  # Set color for normal and big mode
  label = TRUE, lab.pos = "out", lab.size = label_size,  # Move annotations to the top of bars
) + labs(
    y = "Runtime (s)", x = ""
  ) + theme_pubr() + scale_y_sqrt(
    labels = scales::comma, breaks = c(1000, 10000, 100000, 200000, 300000)
  ) + theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = font_size), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  ) + coord_cartesian(clip = 'off')



# Create the ggplot with a secondary y-axis for runtime

fig_1f <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "gray", color = NA))

fig_1g_serine_query_accuracy <- ggbarplot(
    serine_data_long, x = "metric", y = "value",
    fill = "type", palette = tool_colors,
    label = TRUE, lab.pos = "out", lab.size = label_size,
    position = position_dodge(0.8),
    color = NA,  # Removes the bar borders
    ylab = "Metric Value", xlab = "",
    legend = "none" # No legend
    # Re-label the x-axis: Folddisco, pyScoMotif, PDB
) + theme(
    text = element_text(size = font_size),
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off')
 
fig_1h_serine_query_runtime <- ggbarplot(
    serine_data[!is.na(serine_data$runtime),], x = "type", y = "runtime",
    fill = "type", palette = tool_colors,
    label = TRUE, lab.pos = "out", lab.size = label_size,
    color = NA,  # Removes bar borders
    ylab = "Runtime (s)", xlab = "",
    legend = "none" # No legend
) + theme(
    text = element_text(size = font_size),
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off')

# Create the metrics plot for query_len == 3 (Precision, Recall, F1-score)
fig_1i <-  ggplot() + theme_void() + theme(plot.background = element_rect(fill = "gray", color = NA))
fig_1j_zinc_query_3_accuracy <- ggbarplot(
  zinc_data_3_long, x = "metric", y = "value", 
  fill = "type", palette = tool_colors,
  label = TRUE, lab.pos = "out", lab.size = label_size,
  position = position_dodge(0.8),
  color = NA,  # Removes the bar borders
  ylab = "Metric Value", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off')

# Create runtime plot for query_len == 3
fig_1k_zinc_query_3_runtime <- ggbarplot(
  zinc_data_3[!is.na(zinc_data_3$runtime),], x = "type", y = "runtime",
  fill = "type", palette = tool_colors,
  label = TRUE, lab.pos = "out", lab.size = label_size,
  color = NA,  # Removes bar borders
  ylab = "Runtime (s)", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off')


# Create the metrics plot for query_len == 4 (Precision, Recall, F1-score)
fig_1l <-  ggplot() + theme_void() + theme(plot.background = element_rect(fill = "gray", color = NA))
fig_1m_zinc_query_4_accuracy <- ggbarplot(
  zinc_data_4_long, x = "metric", y = "value", 
  fill = "type", palette = tool_colors,
  label = TRUE, lab.pos = "out", lab.size = label_size,
  position = position_dodge(0.8),
  color = NA,  # Removes the bar borders
  ylab = "Metric Value", xlab = "",
  legend = "bottom" # No legend
) + theme(
  text = element_text(size = font_size),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off')

# Create runtime plot for query_len == 4
fig_1n_zinc_query_4_runtime <- ggbarplot(
  zinc_data_4[!is.na(zinc_data_4$runtime),], x = "type", y = "runtime",
  fill = "type", palette = tool_colors,
  label = TRUE, lab.pos = "out", lab.size = label_size,
  color = NA,  # Removes bar borders
  ylab = "Runtime (s)", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off')

# 03. Combine the plots

# Figure 1a will be integrated later
# Just put placeholders for now. fill gray box
fig_1a <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "gray", color = NA))

# Combine the plots into three panels
fig_1_2nd_row <- fig_1b_num_structures + fig_1c_index_size + fig_1d_foldcomp + fig_1e_runtime + 
    plot_layout(ncol = 4, widths = c(1, 2, 1, 2))

# Tagging
# fig_1_2nd_row <- fig_1_2nd_row + plot_annotation(tag_levels = 'a') & theme(
#   plot.tag.position = c(0, 1), plot.tag = element_text(size = 11, face = "bold")
# )


fig_1_3rd_row <- fig_1f + fig_1g_serine_query_accuracy + fig_1h_serine_query_runtime + plot_layout(ncol = 3, widths = c(2, 3.5, 1))
fig_1_4th_row <- fig_1i + fig_1j_zinc_query_3_accuracy + fig_1k_zinc_query_3_runtime + plot_layout(ncol = 3, widths = c(1.5, 3.5, 1))
fig_1_5th_row <- fig_1l + fig_1m_zinc_query_4_accuracy + fig_1n_zinc_query_4_runtime + plot_layout(ncol = 3, widths = c(1.5, 3.5, 1))

fig_1 <- fig_1a / fig_1_2nd_row / fig_1_3rd_row / fig_1_4th_row / fig_1_5th_row + 
  plot_layout(ncol = 1, heights = c(50, 40, 25, 25, 25))


# Testing
ggsave("result/img/fig_1.pdf", fig_1, width = 180, height = 170, unit = "mm", bg = "transparent")







################################################################################


## 99. Supplementary Figures
# Filter the data for h_sapiens, swissprot, pdb, and afdb_cluster_rep with multiple threads
line_plot_data <- index_benchmark_data %>%
  filter(data %in% c('h_sapiens', 'swissprot', 'pdb', 'afdb_cluster_rep')) %>%
  filter (description %in% c('default', 'num_threads', 'database'))

# Create the line plot for num_threads vs runtime with gray color scales
line_plot <- ggplot(line_plot_data, aes(x = num_threads, y = runtime_in_sec, group = data, color = data)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_text(aes(label = round(runtime_in_sec, 2)), vjust = -0.5, size = 4) +
  scale_color_grey(start = 0.2, end = 0.8) +  # Apply gray color scale
  geom_text(
    aes(label = data),
    data = line_plot_data %>% group_by(data) %>% filter(num_threads == max(num_threads)), 
    hjust = -0.2, vjust = 0.5, size = 5) + 
  labs(
    x = "Number of Threads",
    y = "Runtime (s)",
  ) +
  scale_x_continuous(breaks = c(1, 4, 8, 12, 16, 32, 64, 128)) +
  # Make y-axis breaks more dense
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme_pubr() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = 18), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  )

# Print the line plot
print(line_plot)
# Save the plot
ggsave("result/img/num_threads_indexing.pdf", line_plot, width = 9, height = 4, bg = "transparent")
