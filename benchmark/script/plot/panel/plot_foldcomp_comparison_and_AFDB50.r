# Load necessary libraries
library(ggpubr)
library(dplyr)
library(tidyr)
# multi panel plot with patchwork
library(patchwork)

# Folddisco: h_sapiens structure - 0.235GB; index 2.45GB; total 2.685GB
# pyScoMotif: h_sapiens structure - 2.1GB; index 6.1GB; total 8.2GB

data <- data.frame(
  data = c('h_sapiens', 'h_sapiens','h_sapiens', 'h_sapiens'),
  size_in_gb = c(0.235, 2.1, 2.45, 6.1),
  type = c('structure', 'structure', 'index', 'index'),
  tool = c('Folddisco','pyScoMotif', 'Folddisco', 'pyScoMotif')
)
# Calculate the total index size in GB
data <- data %>%
  group_by(data) %>%
  mutate(total_size_in_gb = sum(size_in_gb))

data.tool <- factor(data$tool, levels = c('Folddisco', 'pyScoMotif'))

# Set fill colors for each tool
tool_colors <- c('Folddisco' = "#FFA200", 'pyScoMotif' = "#D3D3D3")

# Create the plot as stacked bar plot
foldcomp_plot <- ggplot(data, aes(x = tool, y = size_in_gb, fill = tool)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  scale_fill_manual(values = tool_colors) +
  labs(y = "Size in GB", x = "") +
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

foldcomp_plot

ggsave("foldcomp_plot.pdf", foldcomp_plot, bg="transparent", width=4.2, height = 4)


# 230	53665860	1379.42

afdb50_data <- data.frame(
  data = c('AFDB50', 'AFDB50', "Swissprot", "Swissprot"),
  size_in_gb = c(230, 1379.42, 3.08, 31.53),
  type = c('structure', 'index', 'structure', 'index')
)

afdb50_data.data <- factor(afdb50_data$data, levels = c('Swissprot', 'AFDB50'))
afdb50_data.colors <- c('Swissprot' = "#D3D3D3", 'AFDB50' = "#FFA200")
  #
afdb50_plot <- ggplot(afdb50_data, aes(x = data, y = size_in_gb, fill = data)) +
  geom_bar(stat = "identity", position = "stack", width = 0.5) +
  scale_fill_manual(values = afdb50_data.colors) +
  labs(y = "Size in GB", x = "") +
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

afdb50_plot
ggsave("afdb50_plot.pdf", afdb50_plot, bg="transparent", width=4.2, height = 4)
