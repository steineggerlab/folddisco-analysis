# Load necessary libraries
library(ggpubr)
library(dplyr)
library(tidyr)
# multi panel plot with patchwork
library(patchwork)


benchmark_data <- read.csv("result/time/index.tsv", header = TRUE, sep = "\t")
# Ensure mode is a factor
benchmark_data$mode <- factor(benchmark_data$index_mode, levels = c("normal", "big"))
benchmark_data$tool <- "FoldDisco"
# Add pyScoMotif and PDB information manually
additional_data <- data.frame(
  data = c('h_sapiens', 'e_coli', 's_cerevisiae', 'pdb', 'pdb'),
  num_structures = c(23391, 4363, 6039, 209564, 209564),
  total_index_size_in_gb = c(6.4, 1.7, 2.2, 73, 55),
  runtime_in_sec = c(3840, 375.89, 775.77, 73800, 298800),  # Converted from HH:MM:SS to seconds
  tool = c('pyScoMotif', 'pyScoMotif', 'pyScoMotif', 'pyScoMotif', 'PDB')
)

# Set fill colors for each tool
tool_colors <- c('FoldDisco' = "#FFA200", 'pyScoMotif' = "#D3D3D3", 'PDB' = "#BEBEBE")



# 1st plot
# Filter for datasets with 12 threads
filtered_data <- benchmark_data %>%
  filter(description %in% c('default', 'other model', 'database') & num_threads == 12) %>%
  arrange(num_structures)
# Combine the data
filtered_data <- bind_rows(filtered_data, additional_data)

# Ensure tool is a factor
filtered_data$tool <- factor(filtered_data$tool, levels = c('FoldDisco', 'pyScoMotif', 'PDB'))



# Create a mapping for the custom x-axis labels
filtered_data$data <- factor(filtered_data$data, 
                             levels = c('e_coli', 's_cerevisiae', 'h_sapiens', 'pdb', 'swissprot'),
                             labels = c('E. coli', 'Yeast', 'Human', 'PDB', 'Swissprot'))

# Create the plot for number of structures
num_structures_plot <- ggbarplot(
  filtered_data %>% distinct(data, num_structures), x = "data", y = "num_structures", fill = "#BEBEBE",
#   palette = c("#323232", "#BEBEBE"),  # Set color for normal and big mode
  label = TRUE, lab.pos = "out", lab.size = 5,  # Move annotations to the top of bars
  color = NA  # Removes bar borders
) +
  scale_y_continuous(breaks = c(100000, 200000, 300000, 400000, 500000, 600000), labels = c("100k", "200k", "300k", "400k", "500k", "600k")) +
  labs(y = "", x = "") +
#   coord_flip() +  # Flip the coordinates
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = 18), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  )

# Create the plot for Index Size
index_size_plot <- ggbarplot(
  filtered_data, x = "data", y = "total_index_size_in_gb", fill = "tool", 
  palette = tool_colors,  # Set color for normal and big mode
  label = TRUE, lab.pos = "out", lab.size = 5,  # Move annotations to the top of bars
  position = position_dodge(),  # Dodge the bars
  color = NA  # Removes bar borders
) +
  labs(
    y = "Index Size (GB)", x = ""
  ) +
  theme_pubr() +
  theme(
     axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = 18), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  )

# Create the plot for Runtime
runtime_plot <- ggbarplot(
  filtered_data, x = "data", y = "runtime_in_sec", fill = "tool", 
  palette = tool_colors,  # Set color for normal and big mode
  label = TRUE, lab.pos = "out", lab.size = 5,  # Move annotations to the top of bars
  position = position_dodge(),  # Dodge the bars
  color = NA  # Removes bar borders
) +
  labs(
    y = "Runtime (s)", x = ""
  ) +
  theme_pubr() +
  scale_y_sqrt(labels = scales::comma, breaks = c(1000, 10000, 100000, 200000, 300000)) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = 18), # Set all font size to 24
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
  )

# Combine the plots into three panels
combined_plot <- num_structures_plot + index_size_plot + runtime_plot + plot_layout(ncol = 3,widths = c(1.8, 3, 3))

# Print the plot
print(combined_plot)

# Save the plot
ggsave("result/img/benchmark_comparison_tools_newer.pdf", combined_plot, width = 18, height = 4, bg = "transparent")

### 1st plot DONE ###
# Filter the data for h_sapiens, swissprot, pdb, and afdb_cluster_rep with multiple threads
line_plot_data <- benchmark_data %>%
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
