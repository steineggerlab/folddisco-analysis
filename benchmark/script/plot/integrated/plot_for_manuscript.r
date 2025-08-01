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

index_benchmark_data <- read.csv("result/indexing/index_pyscomotif_folddisco_t64.tsv", header = TRUE, sep = "\t")
index_benchmark_data_parsed_wopdb <- index_benchmark_data_parsed[!index_benchmark_data_parsed$data %in% c("PDB"),]

index_runtime_data <- read.csv("result/indexing/index_pyscomotif_folddisco_t64.tsv", header = TRUE, sep = "\t")
index_runtime_data$tool <- factor(index_runtime_data$tool, levels = c('Folddisco', 'pyScoMotif', 'pyScoMotif_ext'))


# foldcomp data for Figure 1d
foldcomp_data <- read.csv("result/indexing/foldcomp_comparison.tsv", header = TRUE, sep = "\t")
# Calculate the total index size in GB
# foldcomp_data <- foldcomp_data %>%
#   group_by(data) %>%
#   mutate(total_size_in_gb = sum(size_in_gb))
foldcomp_data$tool <- factor(foldcomp_data$tool, levels = c('Folddisco', 'pyScoMotif', 'MASTER'))

# SCOPe result
pyscomotif_fp1_data <- read.csv("result/SCOPe/250124_scope40_pyscomotif_benchmark_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)
folddisco_prefilter_fp1_data <- read.csv("result/SCOPe/scope40_folddisco_benchmark_result_idf.fp1.parsed.pyscomotif_matched.tsv", sep = "\t", header = TRUE)
folddisco_fp1_data <- read.csv("result/SCOPe/scope40_folddisco_benchmark_result_rmsd.fp1.parsed.pyscomotif_matched.tsv", sep = "\t", header = TRUE)
# folddisco_fp1_data <- read.csv("result/SCOPe/scope40_folddisco_benchmark_result_fixed.tsv_fp1.tsv.parsed.tsv", sep = "\t", header = TRUE)

# Join the two dataframes
pyscomotif_fp1_data$method <- "pyScoMotif"
folddisco_fp1_data$method <- "Folddisco_rmsd"
folddisco_prefilter_fp1_data$method <- "Folddisco_idf"
scope_data <- rbind(pyscomotif_fp1_data,  folddisco_prefilter_fp1_data, folddisco_fp1_data)

scope_colors <- c(
  "pyScoMotif" = "#655AD0", 
  # Different color than Folddisco coloring for RMSD
  "Folddisco_rmsd" = "#FF7F00",
  "Folddisco_idf" = "#E94B8B"
)

# Querying time
querying_time_data <- read.csv("result/querying/querying_time.tsv", header = TRUE, sep = "\t")

# Read and process query_time.tsv
new_querying_time_data <- read.csv("result/querying/query_time.tsv", header = FALSE, sep = "\t")
colnames(new_querying_time_data) <- c("query_file", "index_path", "num_threads", "options", "tool", "runtime", "num_results", "output_file")

# Read and process pyscomotif_time.tsv
pyscomotif_time_data <- read.csv("result/querying/pyscomotif_time.tsv", header = FALSE, sep = "\t")
colnames(pyscomotif_time_data) <- c("query_file", "index_path", "motif_residues", "runtime", "num_results", "output_file")

# Process pyscomotif_time_data to match the structure
pyscomotif_time_processed <- pyscomotif_time_data %>%
  mutate(
    database = case_when(
      grepl("e_coli", index_path) ~ "e_coli",
      grepl("s_cerevisiae", index_path) ~ "s_cerevisiae",
      grepl("h_sapiens", index_path) ~ "h_sapiens", 
      grepl("swissprot", index_path) ~ "swissprot",
      grepl("afdb_rep_v4", index_path) ~ "afdb_cluster_rep",
      grepl("afdb50", index_path) ~ "afdb50",
      TRUE ~ "other"
    ),
    num_structures = as.numeric(case_when(
      database == "e_coli" ~ 4363,
      database == "s_cerevisiae" ~ 6039,
      database == "h_sapiens" ~ 23391,
      database == "swissprot" ~ 542378,
      database == "afdb_cluster_rep" ~ 2266735,
      database == "afdb50" ~ 53665860,
      TRUE ~ 1000
    )),
    tool = "pyScoMotif",
    dataset_source = "pyScoMotif_queries",
    num_threads = 12,
    options = "pyScoMotif_default"
  ) %>%
  select(query_file, index_path, num_threads, options, tool, runtime, num_results, database, num_structures, dataset_source)

# Process existing query_time_data
new_querying_time_data <- new_querying_time_data %>%
  mutate(
    database = case_when(
      grepl("e_coli", index_path) ~ "e_coli",
      grepl("s_cerevisiae", index_path) ~ "s_cerevisiae",
      grepl("h_sapiens", index_path) ~ "h_sapiens", 
      grepl("swissprot", index_path) ~ "swissprot",
      grepl("afdb_rep_v4", index_path) ~ "afdb_cluster_rep",
      grepl("afdb50", index_path) ~ "afdb50",
      TRUE ~ "other"
    ),
    num_structures = as.numeric(case_when(
      database == "e_coli" ~ 4363,
      database == "s_cerevisiae" ~ 6039,
      database == "h_sapiens" ~ 23391,
      database == "swissprot" ~ 542378,
      database == "afdb_cluster_rep" ~ 2266735,
      database == "afdb50" ~ 53665860,
      TRUE ~ 1000
    )),
    dataset_source = "M-CSA_queries"
  ) %>%
  # Ensure tool is a factor with proper levels
  mutate(tool = factor(tool, levels = c('Folddisco', 'pyScoMotif', 'pyScoMotif_ext', 'Folddisco_prefilter'))) %>%
  select(query_file, index_path, num_threads, options, tool, runtime, num_results, database, num_structures, dataset_source)

# Combine both datasets
combined_querying_data <- bind_rows(
  new_querying_time_data,
  pyscomotif_time_processed
) %>%
  mutate(
    tool = factor(tool, levels = c('Folddisco', 'pyScoMotif', 'pyScoMotif_ext', 'Folddisco_prefilter')),
    num_structures = as.numeric(num_structures)
  ) %>%
  group_by(database, tool, num_structures)
# Ensure

# new_querying_time_data$num_structures <- as.numeric(new_querying_time_data$num_structures)
# new_querying_time_data

# query benchmark for Figure 1e - 1l? 
query_benchmark_data <- read.csv("result/querying/query_folddisco.tsv", header = TRUE, sep = "\t")
# Pivot data for easier plotting, excluding accuracy
query_benchmark_data_long <- query_benchmark_data %>%
  select(-accuracy) %>%
  pivot_longer(cols = precision:f1_score, names_to = "metric", values_to = "value")

# Split data for each figure
serine_data <- query_benchmark_data %>% filter(rank == "S01")
serine_data$type <- factor(serine_data$type, levels = c("Folddisco","Folddisco_prefilter", "pyScoMotif", "RCSB", "MASTER"))
zinc_data_3 <- query_benchmark_data %>% filter(query_len == 3 & rank == "C2H2")
zinc_data_3$type <- factor(zinc_data_3$type, levels = c("Folddisco", "Folddisco_prefilter","pyScoMotif", "RCSB", "MASTER"))
zinc_data_4 <- query_benchmark_data %>% filter(query_len == 4 & rank == "C2H2")
zinc_data_4$type <- factor(zinc_data_4$type, levels = c("Folddisco","Folddisco_prefilter", "pyScoMotif", "RCSB", "MASTER"))
serine_data_long <- query_benchmark_data_long %>% filter(rank == "S01")
serine_data_long$type <- factor(serine_data_long$type, levels = c("Folddisco","Folddisco_prefilter", "pyScoMotif", "RCSB", "MASTER"))
zinc_data_3_long <- query_benchmark_data_long %>% filter(query_len == 3 & rank == "C2H2")
zinc_data_3_long$type <- factor(zinc_data_3_long$type, levels = c("Folddisco", "Folddisco_prefilter","pyScoMotif", "RCSB", "MASTER"))
zinc_data_4_long <- query_benchmark_data_long %>% filter(query_len == 4 & rank == "C2H2")
zinc_data_4_long$type <- factor(zinc_data_4_long$type, levels = c("Folddisco", "Folddisco_prefilter","pyScoMotif", "RCSB", "MASTER"))

master_benchmark_data <- read.csv("result/querying/zinc_finger.vs_master.tsv", header = TRUE, sep = "\t")
master_benchmark_data$type <- factor(master_benchmark_data$type, levels = c("Folddisco","Folddisco_prefilter", "pyScoMotif", "RCSB", "MASTER"))
master_benchmark_data_long <- master_benchmark_data %>%
  select(-accuracy) %>%
  pivot_longer(cols = precision:f1_score, names_to = "metric", values_to = "value")
master_benchmark_data_long$type <- factor(master_benchmark_data_long$type, levels = c("Folddisco","Folddisco_prefilter", "pyScoMotif", "RCSB", "MASTER"))

_calc_fit_equation <- function(data, x_col, y_col, tool_name) {
  tool_data <- data[data$tool == tool_name, ]
  if(nrow(tool_data) > 1) {
    model <- lm(get(y_col) ~ get(x_col), data = tool_data)
    intercept <- round(coef(model)[1], 2) # Round to 2 decimal places for intercept
    # For format, use Xe-Y format
    # slope <- round(coef(model)[2], 6)  
    slope <- format(coef(model)[2], scientific = TRUE, digits = 2) # Format slope in scientific notation
    # r_sq <- round(summary(model)$r.squared, 2)
    return(paste0("y = ", slope, "x + ", intercept, "\n"))
  }
  return("")
}

# Fit with exponetial equation
# Power of ?
calc_fit_equation <- function(data, x_col, y_col, tool_name) {
  tool_data <- data[data$tool == tool_name, ]
  if(nrow(tool_data) > 1) {
    # Fit with exponential equation`
    model <- lm(log(get(y_col)) ~ log(get(x_col)), data = tool_data)
    # model <- lm(log(get(y_col)) ~ log(get(x_col)), data = tool_data)
    intercept <- format(exp(coef(model)[1]), digits=2, scientific = TRUE) # Round to 2 decimal places for intercept
    slope <- format(coef(model)[2], digits=2, scientific = TRUE)  # Round to 2 decimal places for slope
    # r_sq <- round(summary(model)$r.squared, 2)
    return(paste0("y = ", intercept, "x^", slope))
  }
  return("")
}


# folddisco_index_size_eq <- calc_fit_equation(index_benchmark_data_parsed_wopdb, "num_structures", "total_index_size_in_gb", "Folddisco")
# folddisco_index_time_eq <- calc_fit_equation(index_runtime_data, "num_structures", "runtime_in_sec", "Folddisco")
# folddisco_querying_time_eq <- calc_fit_equation(querying_time_data, "num_structures", "runtime_in_sec", "Folddisco")
# folddisco_prefilter_querying_time_eq <- calc_fit_equation(querying_time_data, "num_structures", "runtime_in_sec", "Folddisco_prefilter")
# pyscomotif_index_size_eq <- calc_fit_equation(index_benchmark_data_parsed_wopdb, "num_structures", "total_index_size_in_gb", "pyScoMotif")
# pyscomotif_index_time_eq <- calc_fit_equation(index_runtime_data, "num_structures", "runtime_in_sec", "pyScoMotif")
# pyscomotif_querying_time_eq <- calc_fit_equation(querying_time_data, "num_structures", "runtime_in_sec", "pyScoMotif")

# folddisco_index_size_eq
# pyscomotif_index_size_eq
# folddisco_index_time_eq
# pyscomotif_index_time_eq
# folddisco_querying_time_eq
# folddisco_prefilter_querying_time_eq
# pyscomotif_querying_time_eq
# 01-2. Set fill colors for each tool
  # apple2 =  c(
  #   'Folddisco' = "#E94B8B", 'pyScoMotif' = "#655AD0", 'RCSB' = "#4FBDF2", 
  #   'pyScoMotif_ext' = "#655AD0", 'Folddisco_prefilter' = "#F7B4A7",
  #   'MASTER' = "#87B3D7" 
  # ),
  # dark_apple =  c(
  #   'Folddisco' = "#AA5066", 'pyScoMotif' = "#4443A0", 'RCSB' = "#6E89B1", 
  #   'pyScoMotif_ext' = "#4443A0", 'Folddisco_prefilter' = "#D69E55",
  #   'MASTER' = "#5F94D4" 
  # ),
  # ora =  c(
  #   'Folddisco' = "#F0426B", 'pyScoMotif' = "#5A4FCF", 'RCSB' = "#06D6A0", 
  #   'pyScoMotif_ext' = "#5A4FCF", 'Folddisco_prefilter' = "#F68EA6",
  #   'MASTER' = "#FFC43D" 
  # ),
  # neon_pink = c(
  #   'Folddisco' = "#F9025B", 'pyScoMotif' = "#888888", 'RCSB' = "#BBBBBB", 
  #   'pyScoMotif_ext' = "#888888", 'Folddisco_prefilter' = "#F69EC5",
  #   'MASTER' = "#444444" 
  # ),
  # yellow_gray = c(
  #   'Folddisco' = "#FFA200", 'pyScoMotif' = "#888888", 'RCSB' = "#BBBBBB",
  #   'pyScoMotif_ext' = "#888888", 'Folddisco_prefilter' = "#FFA20088",
  #   'MASTER' = "#444444"
  # ),
  # manet = c(
  #   'Folddisco' = "#D29C44", 'pyScoMotif' = "#4585B7", 'RCSB' = "#225E92", 
  #   'pyScoMotif_ext' = "#4585B7", 'Folddisco_prefilter' = "#EBC174",
  #   'MASTER' = "#2F2976"
  # ),
  # neon_hongkong = c(
  #   'Folddisco' = "#CB0C59", 'pyScoMotif' = "#22A0B6", 'RCSB' = "#0B4383", 
  #   'pyScoMotif_ext' = "#22A0B6", 'Folddisco_prefilter' = "#EB64A0",
  #   'MASTER' = "#06394C" 
  # ),
  # palette1 = c(
  #   'Folddisco' = "#E81D52", 'pyScoMotif' = "#043D80", 'RCSB' = "#88BBF0", 
  #   'pyScoMotif_ext' = "#043D80", 'Folddisco_prefilter' = "#FD9D12",
  #   'MASTER' = "#02124E"
  # ),
  # palette2 = c(
  #   'Folddisco' = "#CF152D", 'pyScoMotif' = "#0881C6", 'RCSB' = "#203B7E", 
  #   'pyScoMotif_ext' = "#0881C6", 'Folddisco_prefilter' = "#DF4826",
  #   'MASTER' = "#FEA30C"
  # ),
  # palette3 = c(
  #   'Folddisco' = "#B0016B", 'pyScoMotif' = "#88BBF0", 'RCSB' = "#130959", 
  #   'pyScoMotif_ext' = "#88BBF0", 'Folddisco_prefilter' = "#510346",
  #   'MASTER' = "#F59300"
  # ),
  # palette4 = c(
  #   'Folddisco' = "#B41E6D", 'pyScoMotif' = "#FEB600", 'RCSB' = "#3FA0B3", 
  #   'pyScoMotif_ext' = "#FEB600", 'Folddisco_prefilter' = "#CE968A",
  #   'MASTER' = "#0363A8"
  # ),
  # palette5 = c(
  #   'Folddisco' = "#BB4B51", 'pyScoMotif' = "#13717D", 'RCSB' = "#964D82", 
  #   'pyScoMotif_ext' = "#13717D", 'Folddisco_prefilter' = "#E47236",
  #   'MASTER' = "#9D7F41"
  # )
palettes_list <- list(
  apple =  c(
    'Folddisco' = "#E94B8B", 'pyScoMotif' = "#655AD0", 'RCSB' = "#4FBDF2", 
    'pyScoMotif_ext' = "#655AD0", 'Folddisco_prefilter' = "#F7B4A7",
    'MASTER' = "#B04BB6" 
  )
)

tool_linetypes <- c('Folddisco' = "solid", 'pyScoMotif' = "solid", 'RCSB' = "solid", 'pyScoMotif_ext' = "dashed", 'Folddisco_prefilter' = "solid")
tool_shapes <- c('Folddisco' = 20, 'pyScoMotif' = 20, 'RCSB' = 20, 'pyScoMotif_ext' = 1, 'Folddisco_prefilter' = 20)
tool_alpha <- c('Folddisco' = 1, 'pyScoMotif' = 1, 'RCSB' = 1, 'pyScoMotif_ext' = 1, 'Folddisco_prefilter' = 1)
type_alpha <- c('structure' = 0.5, 'index' = 1)
font_size <- 7
label_size <- 1.6
annotate_size <- 2.0
axis_line_weight <- 0.3
ratio_linetype <- c(
  "1" = "solid",
  "0.6" = "dashed", 
  "0.2" = "dotdash"
)


for (pal_name in names(palettes_list)) {

tool_colors <- palettes_list[[pal_name]]

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
  # annotate(
  #   "text", x = 20000, y = 8000, label = folddisco_index_size_eq, 
  #   size = annotate_size, hjust = 0, vjust = 1, color = tool_colors['Folddisco']
  # ) +
  # annotate(
  #   "text", x = 20000, y = 10000, label = pyscomotif_index_size_eq, 
  #   size = annotate_size, hjust = 0, vjust = 1, color = tool_colors['pyScoMotif']
  # )

# Create the plot for Runtime
indexing_time <- ggline(
  index_runtime_data %>% 
  mutate(runtime_in_sec = round(runtime_in_sec, 0)),
  x = "num_structures", y = "runtime_in_sec", color = "tool", shape = "tool", linetype = "tool",
  palette = tool_colors,  # Set color for normal and big mode
) + labs(
    y = "Runtime for Indexing", x = ""
  ) + theme_pubr() +
  scale_x_sqrt(
    breaks = c(20000, 2000000, 20000000, 50000000),
    labels = c("20K", "2M", "20M", "50M")
  ) +
  scale_y_sqrt(
    breaks = c(3600, 86400, 604800, 1728000),
    labels = c("1h", "1d", "7d", "20d"),
    limit = c(0, 1728000)  # Limit to 20 days
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
  )
  # annotate(
  #   "text", x = 20000, y = 2073600, label = folddisco_index_time_eq,
  #   size = annotate_size, hjust = 0, vjust = 1, color = tool_colors['Folddisco']
  # ) +
  # annotate(
  #   "text", x = 20000, y = 2592000, label = pyscomotif_index_time_eq,
  #   size = annotate_size, hjust = 0, vjust = 1, color = tool_colors['pyScoMotif']
  # )

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


new_querying_time <- ggline(
  combined_querying_data, x = "num_structures", y = "runtime", color = "tool", 
  shape = "tool", linetype = "tool", alpha = "tool",
  palette = tool_colors,  # Set color for normal and big mode
) + labs(
    y = "Runtime for Querying (s)", x = ""
  ) + theme_pubr() +
  scale_shape_manual(values = tool_shapes) +
  scale_linetype_manual(values = tool_linetypes) +
  scale_alpha_manual(values = tool_alpha) +
  scale_x_sqrt(
    breaks = c(20000, 500000, 2000000, 50000000),
    labels = c("20K", "500K", "2M", "50M")
  ) +
  scale_y_continuous(
    breaks = c(0.01, 0.1, 1, 10),
    labels = c("0.01", "0.1", "1", "10"),
    trans = "log10"
  ) +
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
new_querying_time

# Create combined runtime summary from both datasets
runtime_summary <- combined_querying_data %>%
  ungroup() %>%
  group_by(database, tool, num_structures, dataset_source) %>%
  summarise(
    mean_runtime = mean(runtime, na.rm = TRUE),
    sd_runtime = sd(runtime, na.rm = TRUE),
    median_runtime = median(runtime, na.rm = TRUE),
    q25_runtime = quantile(runtime, 0.25, na.rm = TRUE),
    q75_runtime = quantile(runtime, 0.75, na.rm = TRUE),
    min_runtime = min(runtime, na.rm = TRUE),
    max_runtime = max(runtime, na.rm = TRUE),
    n_queries = n(),
    .groups = 'drop'
  ) %>%
  filter(!is.na(tool), !is.na(mean_runtime))

# Create overall summary (combining both datasets for each tool/database)
overall_runtime_summary <- combined_querying_data %>%
  ungroup() %>%
  group_by(database, tool, num_structures) %>%
  summarise(
    mean_runtime = mean(runtime, na.rm = TRUE),
    sd_runtime = sd(runtime, na.rm = TRUE),
    median_runtime = median(runtime, na.rm = TRUE),
    q25_runtime = quantile(runtime, 0.25, na.rm = TRUE),
    q75_runtime = quantile(runtime, 0.75, na.rm = TRUE),
    min_runtime = min(runtime, na.rm = TRUE),
    max_runtime = max(runtime, na.rm = TRUE),
    n_queries = n(),
    datasets_combined = paste(unique(dataset_source), collapse = " + "),
    .groups = 'drop'
  ) %>%
  filter(!is.na(tool), !is.na(mean_runtime)) %>%
  mutate(across(where(is.numeric) & !matches("num_structures|n_queries"), ~ round(.x, 4))) %>%
  mutate(num_structures = as.numeric(num_structures))

# Add extrapolated pyScoMotif data including swissprot for continuous line
extrapolated_pyscomotif <- data.frame(
  database = c("swissprot", "afdb_cluster_rep", "afdb50"),
  tool = c("pyScoMotif_ext", "pyScoMotif_ext", "pyScoMotif_ext"),
  num_structures = c(542378, 2266735, 53665860),
  mean_runtime = c(16.7537, 46.3128, 942.8885),
  sd_runtime = c(10.2374, NA, NA),  # Keep original swissprot error, NA for extrapolated
  median_runtime = c(15.367, 41.1856, 824.183),
  q25_runtime = c(11.2105, NA, NA),  # Keep original swissprot IQR, NA for extrapolated
  q75_runtime = c(20.628, NA, NA),   # Keep original swissprot IQR, NA for extrapolated
  min_runtime = c(0.609, 35.0, 800.0),
  max_runtime = c(113.141, 70.0, 1400.0),
  n_queries = c(751, 751, 751),
  datasets_combined = c("pyScoMotif_extrapolated", "pyScoMotif_extrapolated", "pyScoMotif_extrapolated")
)

# Combine with existing data
overall_runtime_summary <- rbind(overall_runtime_summary, extrapolated_pyscomotif)

# Save to file
write.table(
  overall_runtime_summary, 
  file = "result/querying/overall_runtime_summary.tsv", 
  sep = "\t", row.names = FALSE, quote = FALSE
)

# Update the runtime_plot to exclude error bars for NA values
runtime_plot <- ggplot(
  overall_runtime_summary, 
  aes(x = num_structures, y = mean_runtime, color = tool, shape = tool, linetype = tool)
) +
  geom_line(size = 0.5) +
  geom_point(size = 1.5) +
  # # Add error bars only for non-NA values (excludes extrapolated data)
  # geom_errorbar(
  #   data = subset(overall_runtime_summary, !is.na(q25_runtime) & !is.na(q75_runtime)),
  #   aes(ymin = q25_runtime, 
  #       ymax = q75_runtime, 
  #       color = tool),
  #   width = 0.05, 
  #   size = 0.3
  # ) +
  # Add error bars for standard deviation
  geom_errorbar(
    data = subset(overall_runtime_summary, !is.na(sd_runtime)),
    # assert that mean_runtime - sd_runtime is not negative
    aes(ymin = mean_runtime - 1.96*(sd_runtime/sqrt(n_queries)), 
        ymax = mean_runtime + 1.96*(sd_runtime/sqrt(n_queries)),
    # aes(ymin = pmax(mean_runtime - sd_runtime, min_runtime),
    #     ymax = mean_runtime + sd_runtime, 
        color = tool),
    width = 250,
    size = 0.3
  ) +
  labs(
    y = "Runtime (s)", 
    x = "Number of Structures"
  ) +
  scale_x_sqrt(
    breaks = c(20000, 2000000, 20000000, 50000000),
    labels = c("20K", "2M", "20M", "50M")
  ) +
  scale_y_sqrt(
    breaks = c(1, 10, 100, 1000),
    labels = c("1", "10", "100", "1000")
  ) +
  scale_color_manual(values = tool_colors) +
  scale_shape_manual(values = tool_shapes) +
  scale_linetype_manual(values = tool_linetypes) +
  theme_pubr() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(size = axis_line_weight),
    legend.position = "none",
    plot.title = element_blank(),
    text = element_text(size = font_size),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.box.background = element_rect(fill = "transparent", color = NA)
  )

runtime_plot

  # annotate(
  #   "text", x = 20000, y = 400, label = folddisco_querying_time_eq,
  #   size = annotate_size, hjust = 0, vjust = 1, color = tool_colors['Folddisco']
  # ) +
  # annotate(
  #   "text", x = 20000, y = 300, label = folddisco_prefilter_querying_time_eq,
  #   size = annotate_size, hjust = 0, vjust = 1, color = tool_colors['Folddisco_prefilter']
  # ) +
  # annotate(
  #   "text", x = 20000, y = 500, label = pyscomotif_querying_time_eq,
  #   size = annotate_size, hjust = 0, vjust = 1, color = tool_colors['pyScoMotif']
  # )

# Figure 1d
# Create the plot as stacked bar plot
# foldcomp_plot<- ggplot(foldcomp_data, aes(x = data, y = size_in_gb, fill = tool, alpha = type)) +
#   geom_bar(stat = "identity", position = "stack", width = 0.5) +
foldcomp_plot <- ggbarplot(foldcomp_data[foldcomp_data$data == "swissprot",], x = "tool", y = "size_in_gb", fill = "tool", alpha = "type", color=NA) +
  scale_fill_manual(values = tool_colors) +
  scale_alpha_manual(values = type_alpha) +
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
    legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
    strip.background = element_blank(), strip.text.x = element_blank() # No facet labels
  ) # Change the x-axis labels
foldcomp_plot<- foldcomp_plot+ scale_x_discrete(labels = c("Ours", "pSM", "MASTER"))

# Create the ggplot with a secondary y-axis for runtime

# scope_recall_fp1_plot <- ggviolin(
#     data[data$ratio %in% c(0.2,0.4,0.6,0.8,1.0), ], x = "ratio", y = "recall", fill = "method", 
#     facet.by = "method",
#     add = "mean",
#     ylab = "Sensitivity at 1st FP",
#     xlab = "Ratio",
#     palette = tool_colors,
#     ggtheme = theme_pubr(),
#     legend = "none",
#     size = 0.3,
#     add.params = list(size = 0.2, shape = 20),
#     font.label = list(size = 1.6),
#     ylim = c(0, 1)) + 
#     # stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2), size = 0.5), vjust = -0.25, hjust=-0.7, color = "black") +
#     stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.5, color = "black") +
#     theme(
#       axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
#       axis.line = element_line(linewidth = axis_line_weight),
#       axis.ticks = element_line(size = axis_line_weight),
#       legend.position = "none",  # Removes the legend
#       plot.title = element_blank(),  # Removes the title
#       text = element_text(size = font_size), # Set all font size to 24
#       panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
#       plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
#       legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
#       legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
#       strip.background = element_blank(), strip.text.x = element_blank() # No facet labels
#     ) 


scope_data$method <- factor(scope_data$method, levels = c("Folddisco_rmsd", "Folddisco_idf", "pyScoMotif"))
scope_data$ratio <- factor(scope_data$ratio, levels = c(1.0, 0.8, 0.6, 0.4, 0.2))

# Reverse the curve by 1 - y transformation
scope_suppl_plot <- ggplot(scope_data[scope_data$ratio %in% c(0.2, 0.6, 1.0), ], 
       aes(recall, y = 1 - ..y.., color = method, linetype = ratio, shape = ratio)) + 
       stat_ecdf(geom = "step", size=0.4) + 
       scale_color_manual(values = scope_colors) +
       scale_linetype_manual(values = ratio_linetype) +
      # Need to mask the y-axis to 0-1
       coord_flip(xlim=c(0,1), ylim=c(0,1)) +
       scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
       scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
       theme_pubr() + 
       # Set the x and y labels
        labs(x = "Sensitivity at 1st FP", 
             y = "Fraction of SCOPe queries") +
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

# Reverse the curve by 1 - y transformation
# Draw only for Folddisco_idf and pyScoMotif
scope_recall_fp1_plot <- ggplot(scope_data[scope_data$ratio %in% c(0.2, 0.6, 1.0) & scope_data$method %in% c("Folddisco_idf", "pyScoMotif"), ], 
       aes(recall, y = 1 - ..y.., color = method, linetype = ratio, shape = ratio)) + 
       stat_ecdf(geom = "step", size=0.4) + 
       scale_color_manual(values = scope_colors) +
       scale_linetype_manual(values = ratio_linetype) +
      # Need to mask the y-axis to 0-1
       coord_flip(xlim=c(0,1), ylim=c(0,1)) +
       scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
       scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
       theme_pubr() + 
       # Set the x and y labels
        labs(x = "Sensitivity at 1st FP", 
             y = "Fraction of SCOPe queries") +
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
scope_recall_fp1_plot


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

 
serine_query_runtime <- ggbarplot(
    serine_data[!is.na(serine_data$runtime),], x = "type", y = "runtime",
    fill = "type", palette = tool_colors,
    # label = TRUE, lab.pos = "out", lab.size = label_size,
    color = NA,  # Removes bar borders
    ylab = "Runtime (s)", xlab = "",
    legend = "none" # No legend
) + theme(
    text = element_text(size = font_size),
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks.y = element_line(size = axis_line_weight),
    axis.ticks.x = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + scale_x_discrete(labels = element_blank()
) + scale_y_continuous(limit = c(0, 6.5)
) + coord_cartesian(clip = 'off')

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
zinc_query_3_runtime  <- ggbarplot(
  zinc_data_3[!is.na(zinc_data_3$runtime),], x = "type", y = "runtime",
  fill = "type", palette = tool_colors,
  # label = TRUE, lab.pos = "out", lab.size = label_size,
  color = NA,  # Removes bar borders
  ylab = "Runtime (s)", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  axis.line = element_line(linewidth = axis_line_weight),
  axis.ticks.y = element_line(size = axis_line_weight),
  axis.ticks.x = element_blank(),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA)  # Transparent legend box
) + scale_x_discrete(labels = element_blank()
) + scale_y_continuous(limit = c(0, 6.5)
) + coord_cartesian(clip = 'off')


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
  legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
) + coord_cartesian(clip = 'off')  + scale_x_discrete(labels = c("Precision", "Recall", "F1-score"))

# # Create runtime plot for query_len == 4
zinc_query_4_runtime <- ggbarplot(
  zinc_data_4[!is.na(zinc_data_4$runtime),], x = "type", y = "runtime",
  fill = "type", palette = tool_colors,
  # label = TRUE, lab.pos = "out", lab.size = label_size,
  color = NA,  # Removes bar borders
  ylab = "Runtime (s)", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  axis.line = element_line(linewidth = axis_line_weight),
  axis.ticks.y = element_line(size = axis_line_weight),
  axis.ticks.x = element_blank(),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
) + scale_x_discrete(
  labels  = element_blank()
) + scale_y_continuous(limit = c(0, 6.5)
) + coord_cartesian(clip = 'off')

master_accuracy <- ggbarplot(
  master_benchmark_data_long, x = "metric", y = "value", 
  fill = "type", palette = tool_colors,
  # label = TRUE, lab.pos = "out", lab.size = label_size,
  position = position_dodge(0.8),
  color = NA,  # Removes the bar borders
  ylab = "Metric Value", xlab = "",
  legend = "none" # No legend
) + scale_y_continuous(limits = c(0, 1)) +
theme(
  text = element_text(size = font_size),
  axis.line = element_line(linewidth = axis_line_weight),
  axis.ticks = element_line(size = axis_line_weight),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
) + coord_cartesian(clip = 'off') + scale_x_discrete(labels = c("Precision", "Recall", "F1-score"))

master_runtime <- ggbarplot(
  master_benchmark_data[!is.na(master_benchmark_data$runtime),], x = "type", y = "runtime",
  fill = "type", palette = tool_colors,
  # label = TRUE, lab.pos = "out", lab.size = label_size,
  color = NA,  # Removes bar borders
  ylab = "Runtime (s)", xlab = "",
  legend = "none" # No legend
) + theme(
  text = element_text(size = font_size),
  axis.line = element_line(linewidth = axis_line_weight),
  axis.ticks.y = element_line(size = axis_line_weight),
  axis.ticks.x = element_blank(),
  panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
  plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
  legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
  legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
) + scale_x_discrete(
  labels  = element_blank()
) + scale_y_continuous(limit = c(0, 300)
) + coord_cartesian(clip = 'off')

stacked_plot <- indexing_time / index_size / querying_time 

# 03. Combine the plots

# Figure 1a will be integrated later
# Just put placeholders for now. fill gray box
#fig_1a <- ggplot() + theme_void() + theme(plot.background = element_rect(fill = "gray", color = NA))

# Combine the plots into three panels

# layout_design <- "############
#                   AA#BBBCCCDDD
#                   EEE###FFFFFF
#                   ##GGGH##IIIJ"
# layout_design <- "############
#                   AA#BBBCCCDDD
#                   ##EEEF##GGGH
#                   III###JJJJJJ"
layout_design <- "##############
                  AAABBBCCCDDDDD
                  EEEE##FFFF####"

layout_design <- "##############
                  AAAA##BBBB####
                  CCCDDDEEEFFFFF"


## Current panels
# 1st row: 1a:num_structures_table, 1b:index_size, 1c:indexing_time, 1d: query_time
# 2nd row: 1e:foldcomp_plot, 1f:scope_recall_fp1_plot, 1g:esm_barplot
# 3rd row: 1h: serine_motif (empty panel), 1i: serine_query_accuracy, 1j: serine_query_runtime, 
# 3rd row (continued) 1k: zinc3 motif (empty panel), 1l: zinc3_query_accuracy, 1m: zinc3_query_runtime
# 4th row: 1n: zinc4 motif (empty panel), 1o: zinc4_query_accuracy, 1p: zinc4_query_runtime
# 4th row (continued) 1q: zinc_patch (empty panel), 1r: zinc_patch_query_accuracy, 1s: zinc_patch_query_runtime 
# Patch query benchmark contains comparison against MASTER

## Updated panels
# 1st row: 1b:index_size, 1c:indexing_time, 1d: query_time
# 2nd row: 1e:foldcomp_plot, 1f:scope_recall_fp1_plot,
# 3rd row: 1n: zinc4 motif (empty panel), 1o: zinc4_query_accuracy, 1p: zinc4_query_runtime
# 3rd row (continued) 1q: zinc_patch (empty panel), 1r: zinc_patch_query_accuracy, 1s: zinc_patch_query_runtime 
# Patch query benchmark contains comparison against MASTER


# TODO: need alignment
# fig2 <- num_structures_table + index_size + indexing_time + querying_time +
#         zinc_query_4_accuracy + zinc_query_4_runtime + master_accuracy + master_runtime +
#         foldcomp_plot + scope_recall_fp1_plot +
#          plot_layout(
#           design = layout_design, 
#           widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#           heights = c(1, 1, 1, 1)
#         )

# fig2 <- index_size + indexing_time + querying_time + scope_recall_fp1_plot +
#         zinc_query_4_accuracy + master_accuracy + 
#         plot_layout(
#           design = layout_design, 
#           widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
#           heights = c(1, 1, 1)
#         )
fig2 <- zinc_query_4_accuracy + master_accuracy + 
        index_size + indexing_time + runtime_plot + scope_recall_fp1_plot +
        plot_layout(
          design = layout_design, 
          widths = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
          heights = c(1, 1, 1)
        )


# Testing
ggsave(sprintf("result/img/fig2_%s.pdf", pal_name), fig2, width = 180, height = 110, unit = "mm", bg = "transparent")

ggsave(sprintf("result/img/fig2_%s.png", pal_name), fig2, width = 180, height = 110, unit = "mm", bg = "transparent")

# Suppl. fig
suppl_layout <- "##AAAB
                 ##CCCD"

suppl_fig1 <- serine_query_accuracy + serine_query_runtime + zinc_query_3_accuracy + zinc_query_3_runtime +
              plot_layout(
                design = suppl_layout,
                widths = c(1, 1, 1, 1),
                heights = c(1, 1)
              )
ggsave(sprintf("result/img/suppl_fig1_%s.pdf", pal_name), suppl_fig1, width = 90, height = 80, unit = "mm", bg = "transparent")
ggsave(sprintf("result/img/suppl_fig1_%s.png", pal_name), suppl_fig1, width = 90, height = 80, unit = "mm", bg = "transparent")



# Suppl fig 2
ggsave("result/img/suppl_fig2.pdf", scope_suppl_plot, width = 90, height = 60, unit = "mm", bg = "transparent")
ggsave("result/img/suppl_fig2.png", scope_suppl_plot, width = 90, height = 60, unit = "mm", bg = "transparent")

# Suppl fig 3




}
################################################################################


# ## 99. Supplementary Figures
# # Filter the data for h_sapiens, swissprot, pdb, and afdb_cluster_rep with multiple threads
line_plot_data <- index_benchmark_data %>%
  filter(data %in% c('e_coli', 'h_sapiens', 'swissprot', 'pdb', 'afdb_cluster_rep', 'afdb50')) %>%
  filter (description %in% c('default', 'num_threads', 'database'))

# # Create the line plot for num_threads vs runtime with gray color scales
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

line_plot

# Create the new data frame with only h_sapiens, swissprot, and e_coli
line_plot_data_filtered <- data.frame(
  data = c("h_sapiens_filtered", "e_coli_filtered", "afdb_swissprot_v4", 
           "h_sapiens_filtered", "e_coli_filtered", "afdb_swissprot_v4", 
           "h_sapiens_filtered", "e_coli_filtered", "afdb_swissprot_v4", 
           "h_sapiens_filtered", "e_coli_filtered", "afdb_swissprot_v4"),
  num_threads = c(64, 64, 64, 16, 16, 16, 4, 4, 4, 1, 1, 1),
  runtime_in_sec = c(52.564, 6.699, 1081.648, 93.846, 8.928, 2381.752, 
                     211.064, 17.663, 5774.601, 728.110, 55.615, 19065.781),
  num_structures = c(23391, 4363, 542378, 23391, 4363, 542378, 
                     23391, 4363, 542378, 23391, 4363, 542378),
)

# Create the filtered line plot
line_plot_filtered <- ggplot(line_plot_data_filtered, aes(x = num_threads, y = runtime_in_sec, group = data, color = data)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_text(aes(label = round(runtime_in_sec, 1)), vjust = -0.5, size = 3) +
  scale_color_grey(start = 0.2, end = 0.8) +  # Apply gray color scale
  geom_text(
    aes(label = data),
    data = line_plot_data_filtered %>% group_by(data) %>% filter(num_threads == max(num_threads)), 
    hjust = -0.1, vjust = 0.5, size = 4
  ) + 
  labs(
    x = "Number of Threads",
    y = "Runtime (s)"
  ) +
  scale_x_continuous(breaks = c(1, 4, 16, 64)) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                # Label format as integers. 1, 10, 100, 1000, etc.
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  theme_pubr() +
  theme(
    legend.title = element_blank(),
    legend.position = "none",  # Removes the legend
    plot.title = element_blank(),  # Removes the title
    text = element_text(size = font_size), 
    panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
    plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
    legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
    legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
    axis.line = element_line(linewidth = axis_line_weight),
    axis.ticks = element_line(size = axis_line_weight)
  )

line_plot_filtered

# # Print the line plot
# print(line_plot)
# # Save the plot
ggsave("result/img/num_threads_indexing.pdf", line_plot_filtered, width = 90, height = 60, units = "mm", bg = "transparent")



# # Create the plot for scope recall at 1st FP
# scope_recall_fp1_plot <- ggecdf(
#     data[data$ratio %in% c(0.2,0.4,0.6,0.8,1.0), ], x = "recall", color = "method", linetype = "ratio", 
#     # facet.by = "method",
#     # add = c("median_q1q3"),

#     # ylab = "Sensitivity at 1st FP",
#     # xlab = "Ratio",
#     palette = tool_colors,
#     ggtheme = theme_pubr(),
#     legend = "none",
#     size = 0.3,
#     add.params = list(size = 0.2, shape = 20),
#     font.label = list(size = 1.6),
#     # Set dotplot parameters
#     ylim = c(0, 1)) + 
#     # stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2), size = 0.5), vjust = -0.25, hjust=-0.7, color = "black") +
#     # stat_summary(fun = mean, geom = "line", aes(group = 1), linewidth = 0.5, color = "black") +
#     theme(
#       axis.text.x = element_text(angle = 0, hjust = 0.5),  # Horizontal axis text
#       axis.line = element_line(linewidth = axis_line_weight),
#       axis.ticks = element_line(size = axis_line_weight),
#       legend.position = "none",  # Removes the legend
#       plot.title = element_blank(),  # Removes the title
#       text = element_text(size = font_size), # Set all font size to 24
#       panel.background = element_rect(fill = "transparent", color = NA),  # Transparent panel
#       plot.background = element_rect(fill = "transparent", color = NA),   # Transparent plot area
#       legend.background = element_rect(fill = "transparent", color = NA), # Transparent legend
#       legend.box.background = element_rect(fill = "transparent", color = NA),  # Transparent legend box
#       strip.background = element_blank(), strip.text.x = element_blank() # No facet labels
#     )

# scope_recall_fp1_plot


# share x-axis. Applying same scale to x-axis to all plots
# Stack
stacked_plot
auc_data <- scope_data[scope_data$ratio %in% c(0.2, 0.6, 1.0), ]

# Function to calculate AUC using trapezoidal rule
calculate_auc <- function(x, y) {
  # Sort by x values
  ord <- order(x)
  x <- x[ord]
  y <- y[ord]
  
  # Remove duplicates and NAs
  valid <- !is.na(x) & !is.na(y)
  x <- x[valid]
  y <- y[valid]
  
  # Calculate AUC using trapezoidal rule
  if(length(x) < 2) return(NA)
  
  # Ensure x goes from 0 to 1 for proper AUC calculation
  if(min(x) > 0) {
    x <- c(0, x)
    y <- c(0, y)
  }
  if(max(x) < 1) {
    x <- c(x, 1)
    y <- c(y, tail(y, 1))
  }
  
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
}

# Calculate AUC for each method and ratio
auc_results <- auc_data %>%
  group_by(method, ratio) %>%
  summarise(
    auc = calculate_auc(recall, ecdf(recall)(recall)),
    .groups = 'drop'
  )

auc_results

auc_results_ecdf <- auc_data %>%
  group_by(method, ratio) %>%
  summarise(
    # AUC = integral of (1 - ECDF(x)) from 0 to 1
    auc = integrate(function(x) 1 - ecdf(recall)(x), 0, 1)$value,
    mean_recall = mean(recall, na.rm = TRUE),
    median_recall = median(recall, na.rm = TRUE),
    .groups = 'drop'
  )

auc_results_ecdf

calculate_auroc1 <- function(data) {
  # For each method and ratio combination
  auroc1_results <- data %>%
    group_by(method, ratio) %>%
    summarise(
      # AUROC1 is the area under the ROC curve up to the first false positive
      # This is equivalent to the sensitivity (recall) at the first FP
      auroc1 = mean(recall, na.rm = TRUE),
      median_recall = median(recall, na.rm = TRUE),
      q25_recall = quantile(recall, 0.25, na.rm = TRUE),
      q75_recall = quantile(recall, 0.75, na.rm = TRUE),
      n_queries = n(),
      .groups = 'drop'
    )
  
  return(auroc1_results)
}

# Calculate AUROC1 for your data
auroc1_scope_results <- calculate_auroc1(scope_data[scope_data$ratio %in% c(0.2, 0.6, 1.0), ])

# Display results
print("AUROC1 Results by Method and Ratio:")
print(auroc1_scope_results)


if (!require(gridExtra)) {
  install.packages("gridExtra")
  library(gridExtra)
}
# Alternative approach using a simpler table creation method
auroc1_table_simple <- ggtexttable(
  auroc1_table_data,
  rows = NULL,
  theme = ttheme("blank", base_size = font_size)
)

# Or try without the formatting first
auroc1_table_basic <- ggtexttable(auroc1_table_data)

# Check if the basic version works
print(auroc1_table_basic)

# If it works, then add formatting step by step
auroc1_table <- auroc1_table_basic %>% 
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)

auroc1_table
auro