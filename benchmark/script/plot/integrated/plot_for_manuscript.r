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
index_benchmark_data_parsed <- index_benchmark_data %>%
  filter(description %in% c('default', 'other_model', 'database') & num_threads == 12) %>%
  arrange(num_structures)
# Combine the data
index_benchmark_data_parsed <- bind_rows(index_benchmark_data_parsed, additional_data)
# Ensure tool is a factor
index_benchmark_data_parsed$tool <- factor(index_benchmark_data_parsed$tool, levels = c('Folddisco', 'pyScoMotif', 'pyScoMotif_ext'))
# Create a mapping for the custom x-axis labels
index_benchmark_data_parsed$data <- factor(index_benchmark_data_parsed$data, 
                             levels = c('e_coli', 'h_sapiens', 'pdb', 'swissprot', 'afdb_cluster_rep', 'afdb50'),
                             labels = c('E. coli', 'Human', 'PDB', 'SwissProt', 'AFDB.FS', 'AFDB50'))

index_benchmark_data_parsed_wopdb <- index_benchmark_data_parsed[!index_benchmark_data_parsed$data %in% c("PDB"),]
index_benchmark_data_parsed_wopdb
# foldcomp data for Figure 1d
foldcomp_data <- read.csv("result/indexing/foldcomp_comparison.tsv", header = TRUE, sep = "\t")
# Calculate the total index size in GB
foldcomp_data <- foldcomp_data %>%
  group_by(data) %>%
  mutate(total_size_in_gb = sum(size_in_gb))
foldcomp_data.tool <- factor(foldcomp_data$tool, levels = c('Folddisco', 'pyScoMotif'))

# SCOPe result
pyscomotif_fp1_data <- read.csv("result/SCOPe/250124_scope40_pyscomotif_benchmark_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)
# folddisco_fp1_data <- read.csv("result/SCOPe/250214_scope_benchmark_with_length_penalty_0.5_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)
folddisco_fp1_data <- read.csv("result/SCOPe/250216_scope_benchmark_with_lp_and_match_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)

# Join the two dataframes
pyscomotif_fp1_data$method <- "pyScoMotif"
folddisco_fp1_data$method <- "Folddisco"
data <- rbind(pyscomotif_fp1_data, folddisco_fp1_data)

# Querying time
querying_time_data <- read.csv("result/querying/querying_time.tsv", header = TRUE, sep = "\t")


# 01-2. Set fill colors for each tool
# Yellow
# tool_colors <- c('Folddisco' = "#FFA200", 'pyScoMotif' = "#AAAAAA", 'RCSB' = "#666666", 'pyScoMotif_ext' = "#AAAAAA", 'Folddisco_skip' = "#FFA20088")
# Neon pink
tool_colors <- c(
  'Folddisco' = "#F9025B", 'pyScoMotif' = "#AAAAAA", 'RCSB' = "#666666", 
  'pyScoMotif_ext' = "#AAAAAA", 'Folddisco_skip' = "#F9025B88",
  'MASTER' = "#222222"
)
# Disco palette from pinterest
# tool_colors <- c('Folddisco' = "#EE227D", 'pyScoMotif' = "#FD8083", 'RCSB' = "#498099", 'pyScoMotif_ext' = "#FD8083", 'Folddisco_skip' = "#EE227D88")
tool_linetypes <- c('Folddisco' = "solid", 'pyScoMotif' = "solid", 'RCSB' = "solid", 'pyScoMotif_ext' = "dashed", 'Folddisco_skip' = "solid")
tool_shapes <- c('Folddisco' = 20, 'pyScoMotif' = 20, 'RCSB' = 20, 'pyScoMotif_ext' = 1, 'Folddisco_skip' = 20)
tool_alpha <- c('Folddisco' = 1, 'pyScoMotif' = 1, 'RCSB' = 1, 'pyScoMotif_ext' = 1, 'Folddisco_skip' = 0.5)
font_size <- 7
label_size <- 1.6
axis_line_weight <- 0.3

# 02. Plotting
# Figure 1a
# Create the plot for number of structures
# ggtexttable
num_structures_table <- ggtexttable(
  index_benchmark_data_parsed %>% distinct(data, num_structures),
  rows = NULL, cols = c("DB", "Structures"),
  theme = ttheme(
    base_style = "blank", base_size = font_size,
    padding = unit(c(0.8, 2), "mm")
)) %>% tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 2)

# Figure 1c
# Create the plot for Index Size
# Currently, using log-log plot. Transform both x and y axis to log10
index_size <- ggline(
  index_benchmark_data_parsed_wopdb, x = "num_structures", y = "total_index_size_in_gb", color = "tool", linetype = "tool", shape = "tool",
  palette = tool_colors,  # Set color for normal and big mode
  # Changed to line; No labels
  label.nb.digits = 1,
) +
  labs(
    y = "Index Size (GB)", x = ""
  ) +
  scale_x_sqrt(
    breaks = c(20000, 2000000, 20000000, 50000000),
    labels = c("20K", "2M", "20M", "50M")
  ) +
  scale_y_sqrt(
    breaks = c(1, 100, 1000, 10000),
    limit = c(0, 10000)
  ) +
  scale_linetype_manual(values = tool_linetypes) + scale_shape_manual(values = tool_shapes) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(size = axis_line_weight),
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = font_size), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
    legend.title = element_blank()
    # Set line weight to axis
  ) 

# Create the plot for Runtime
indexing_time <- ggline(
  index_benchmark_data_parsed_wopdb %>% 
  mutate(runtime_in_sec = round(runtime_in_sec, 0)) %>%
  filter(data != "AFDB50"),
  x = "num_structures", y = "runtime_in_sec", color = "tool", shape = "tool", linetype = "tool",
  palette = tool_colors,  # Set color for normal and big mode
) + labs(
    y = "Runtime for Indexing", x = ""
  ) + theme_pubr() +
  scale_x_sqrt(
    breaks = c(20000, 500000, 2000000),
    labels = c("20K", "500K", "2M")
  ) +
  scale_y_sqrt(
    breaks = c(60, 3600, 28800, 86400, 259200),
    labels = c("1m", "1h", "8h", "1d", "3d"),
    limit = c(0, 300000)
  ) +
  scale_linetype_manual(values = tool_linetypes) + scale_shape_manual(values = tool_shapes) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(size = axis_line_weight),
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = font_size), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  ) 

# Fig 1d --> quey speed
querying_time <- ggline(
  querying_time_data, x = "num_structures", y = "runtime_in_sec", color = "tool", 
  shape = "tool", linetype = "tool", alpha = "tool",
  palette = tool_colors,  # Set color for normal and big mode
) + labs(
    y = "Runtime for Querying (s)", x = ""
  ) + theme_pubr() +
  scale_shape_manual(values = tool_shapes) +
  scale_linetype_manual(values = tool_linetypes) +
  scale_alpha_manual(values = tool_alpha) +
  scale_x_sqrt(
    breaks = c(20000, 2000000, 20000000, 50000000),
    labels = c("20K", "2M", "20M", "50M")
  ) +
  scale_y_sqrt(
    breaks = c(1, 10, 100, 200, 500),
    labels = c("1", "10", "100", "200", "500")
  ) +
  # Not using log scale
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(size = axis_line_weight),
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = font_size), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  )


# Figure 1d
# Create the plot as stacked bar plot
foldcomp_plot<- ggplot(foldcomp_data, aes(x = tool, y = size_in_gb, fill = tool, alpha = type)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  scale_fill_manual(values = tool_colors) +
  scale_alpha_manual(values = c(0.5, 1)) +
  labs(y = "Size in GB", x = "") +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(size = axis_line_weight),
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = font_size), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  ) # Change the x-axis labels
foldcomp_plot<- foldcomp_plot+ scale_x_discrete(labels = c("Ours", "pyScoMotif"))

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

master_benchmark_data <- read.csv("result/querying/zinc_finger.vs_master.tsv", header = TRUE, sep = "\t")
master_benchmark_data_long <- master_benchmark_data %>%
  select(-accuracy) %>%
  pivot_longer(cols = precision:f1_score, names_to = "metric", values_to = "value")
master_benchmark_data_long$type <- factor(master_benchmark_data_long$type, levels = c("Folddisco", "pyScoMotif", "RCSB", "MASTER"))


# Create the ggplot with a secondary y-axis for runtime


scope_recall_fp1_plot <- ggviolin(
    data[data$ratio %in% c(0.2,0.4,0.6,0.8,1.0), ], x = "ratio", y = "recall", fill = "method", 
    facet.by = "method",
    add = "mean",
    ylab = "Sensitivity at 1st FP",
    xlab = "Ratio",
    palette = tool_colors,
    ggtheme = theme_pubr(),
    legend = "none",
    size = 0.3,
    add.params = list(size = 0.2, shape = 20),
    font.label = list(size = 1.6),
    ylim = c(0, 1)) + 
    # stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2), size = 0.5), vjust = -0.25, hjust=-0.7, color = "black") +
    stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.5, color = "black") +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
      axis.line = element_line(linewidth = axis_line_weight),
      axis.ticks = element_line(size = axis_line_weight),
      legend.position = "none",  # Removes the legend
      plot.title = element_blank(),  # Removes the title
      text = element_text(size = font_size), # Set all font size to 24
      panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
      plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
      legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
      legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
      strip.background = element_blank(), strip.text.x = element_blank() # No facet labels
    ) 


serine_query_accuracy <- ggbarplot(
    serine_data_long, x = "metric", y = "value",
    fill = "type", palette = tool_colors,
    # label = TRUE, lab.pos = "out", lab.size = label_size,
    position = position_dodge(0.8),
    color = NA,  # Removes the bar borders
    ylab = "Metric Value", xlab = "",
    legend = "none" # No legend
    # Re-label the x-axis: Folddisco, pyScoMotif, PDB
) + theme(
    text = element_text(size = font_size),
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(size = axis_line_weight),
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off') + scale_x_discrete(labels = c("Precision", "Recall", "F1-score"))

 
# serine_query_runtime <- ggbarplot(
#     serine_data[!is.na(serine_data$runtime),], x = "type", y = "runtime",
#     fill = "type", palette = tool_colors,
#     # label = TRUE, lab.pos = "out", lab.size = label_size,
#     color = NA,  # Removes bar borders
#     ylab = "Runtime (s)", xlab = "",
#     legend = "none" # No legend
# ) + theme(
#     text = element_text(size = font_size),
#     axis.line = element_line(linewidth = axis_line_weight),
#     axis.ticks.y = element_line(size = axis_line_weight),
#     axis.ticks.x = element_blank(),
#     panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
#     plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
#     legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
#     legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
# ) + scale_x_discrete(labels = element_blank()
# ) + scale_y_continuous(limit = c(0, 6.5)
# ) + coord_cartesian(clip = 'off')
# serine_query_runtime
# Create the metrics plot for query_len == 3 (Precision, Recall, F1-score)
zinc_query_3_accuracy  <- ggbarplot(
  zinc_data_3_long, x = "metric", y = "value", 
  fill = "type", palette = tool_colors,
  # label = TRUE, lab.pos = "out", lab.size = label_size,
  position = position_dodge(0.8),
  color = NA,  # Removes the bar borders
  ylab = "Metric Value", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  axis.line = element_line(linewidth = axis_line_weight),
  axis.ticks = element_line(size = axis_line_weight),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off') + scale_x_discrete(labels = c("Precision", "Recall", "F1-score"))

# # Create runtime plot for query_len == 3
# zinc_query_3_runtime  <- ggbarplot(
#   zinc_data_3[!is.na(zinc_data_3$runtime),], x = "type", y = "runtime",
#   fill = "type", palette = tool_colors,
#   # label = TRUE, lab.pos = "out", lab.size = label_size,
#   color = NA,  # Removes bar borders
#   ylab = "Runtime (s)", xlab = "",
#   legend = "none" # No legend
# ) + theme(
#   text = element_text(size = font_size),
#   axis.line = element_line(linewidth = axis_line_weight),
#   axis.ticks = element_line(size = axis_line_weight),
#   panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
#   plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
#   legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
#   legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
# ) + scale_x_discrete(labels = element_blank()
# ) + scale_y_continuous(limit = c(0, 6.5)
# ) + coord_cartesian(clip = 'off')


# Create the metrics plot for query_len == 4 (Precision, Recall, F1-score)
zinc_query_4_accuracy <- ggbarplot(
  zinc_data_4_long, x = "metric", y = "value", 
  fill = "type", palette = tool_colors,
  # label = TRUE, lab.pos = "out", lab.size = label_size,
  position = position_dodge(0.8),
  color = NA,  # Removes the bar borders
  ylab = "Metric Value", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  axis.line = element_line(linewidth = axis_line_weight),
  axis.ticks = element_line(size = axis_line_weight),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off')  + scale_x_discrete(labels = c("Precision", "Recall", "F1-score"))

# # Create runtime plot for query_len == 4
# zinc_query_4_runtime <- ggbarplot(
#   zinc_data_4[!is.na(zinc_data_4$runtime),], x = "type", y = "runtime",
#   fill = "type", palette = tool_colors,
#   # label = TRUE, lab.pos = "out", lab.size = label_size,
#   color = NA,  # Removes bar borders
#   ylab = "Runtime (s)", xlab = "",
#   legend = "none" # No legend
# ) + theme(
#   text = element_text(size = font_size),
#   axis.line = element_line(linewidth = axis_line_weight),
#   axis.ticks = element_line(size = axis_line_weight),
#   panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
#   plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
#   legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
#   legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
# ) + scale_x_discrete(
#   labels  = element_blank()
# ) + scale_y_continuous(limit = c(0, 6.5)
# ) + coord_cartesian(clip = 'off')

master_accuracy <- ggbarplot(
  master_benchmark_data_long, x = "metric", y = "value", 
  fill = "type", palette = tool_colors,
  # label = TRUE, lab.pos = "out", lab.size = label_size,
  position = position_dodge(0.8),
  color = NA,  # Removes the bar borders
  ylab = "Metric Value", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  axis.line = element_line(linewidth = axis_line_weight),
  axis.ticks = element_line(size = axis_line_weight),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + coord_cartesian(clip = 'off') + scale_x_discrete(labels = c("Precision", "Recall", "F1-score"))

# 03. Combine the plots

# Figure 1a will be integrated later
# Just put placeholders for now. fill gray box
#fig_1a <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "gray", color = NA))

# Combine the plots into three panels

layout_design <- "AA#BBBCCCDDD
                  EEEFFFFFF###
                  ##GGGG##HHHH
                  ##IIII##JJJJ"

## Current panels
# 1st row: 1a:num_structures_table, 1b:index_size, 1c:indexing_time, 1d: query_time
# 2nd row: 1e:foldcomp_plot, 1f:scope_recall_fp1_plot, 1g:esm_barplot
# 3rd row: 1h: serine_motif (empty panel), 1i: serine_query_accuracy, 1j: serine_query_runtime, 
# 3rd row (continued) 1k: zinc3 motif (empty panel), 1l: zinc3_query_accuracy, 1m: zinc3_query_runtime
# 4th row: 1n: zinc4 motif (empty panel), 1o: zinc4_query_accuracy, 1p: zinc4_query_runtime
# 4th row (continued) 1q: zinc_patch (empty panel), 1r: zinc_patch_query_accuracy, 1s: zinc_patch_query_runtime 
# Patch query benchmark contains comparison against MASTER

# TODO: need alignment
fig2 <- num_structures_table + index_size + indexing_time + querying_time +
         foldcomp_plot+ scope_recall_fp1_plot +
         serine_query_accuracy + zinc_query_3_accuracy +
         zinc_query_4_accuracy + master_accuracy +
         plot_layout(
          design = layout_design, 
          widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
          heights = c(1, 1, 1, 1, 1)
        )

# Testing
ggsave("result/img/fig2.pdf", fig2, width = 180, height = 160, unit = "mm", bg = "transparent")





################################################################################


## 99. Supplementary Figures
# Filter the data for h_sapiens, swissprot, pdb, and afdb_cluster_rep with multiple threads
line_plot_data <- index_benchmark_data %>%
  filter(data %in% c('e_coli', 'h_sapiens', 'swissprot', 'pdb', 'afdb_cluster_rep', 'afdb50')) %>%
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
