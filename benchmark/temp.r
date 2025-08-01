# Print heatmap
library(pheatmap)
library(ggpubr)
# Plasma palette
library(viridis)


data <- read.csv("benchmark/AF-O15391-F1-model_v4.pdb.s_cerevisiae.idf_matrix.csv", header = FALSE)
data <- read.csv("benchmark/AF-O15391-F1-model_v4.pdb.h_sapiens.idf_matrix.csv", header = FALSE)
data <- read.csv("benchmark/AF-O15391-F1-model_v4.h_sapiens.hybrid.idf_matrix.csv", header = FALSE)

data <- read.csv("/Users/hbk/Projects/Lab/06_FoldMotif/repos/folddisco_analysis_organizing/folddisco-analysis/benchmark/result/idf_testing/serine_protease_ppf.idf_matrix.csv", header = FALSE)

data <- read.csv("/Users/hbk/Projects/Lab/06_FoldMotif/repos/folddisco_analysis_organizing/folddisco-analysis/benchmark/result/idf_testing/zinc_finger_folddisco.idf_matrix.csv", header = FALSE)


data <- as.matrix(data)

subdata <- data[c(20:50),c(20:50)]

subdata <- data[c(10:40),c(60:100)]
subdata
pheatmap(data, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = FALSE, show_colnames = FALSE, color = magma(16))


vec_from_data <- function(data) {
  vec <- c()
  for (i in 1:dim(data)[1]) {
    for (j in 1:dim(data)[2]) {
        # Skip 0
        if (data[i, j] == 0) {
            next
        }
        vec <- c(vec, data[i, j])
    }
  }
  return(vec)
}

vec_data <- vec_from_data(data)
vec_sub_data <- vec_from_data(subdata)
motif_data <- read.csv("/Users/hbk/Projects/Lab/06_FoldMotif/repos/folddisco_analysis_organizing/folddisco-analysis/benchmark/result/idf_testing/zinc_finger_folddisco.motif_idf_vector.csv", header = FALSE)

# Overlay density plot
vec_motif_data <- vec_from_data(as.matrix(motif_data))

df <- data.frame(
  idf_value = c(vec_data, vec_sub_data),
  group = c(
    rep("Data", length(vec_data)),
    rep("Helix", length(vec_sub_data))
  )
)

ggdensity(df, x = "idf_value", fill = "group", palette = c("blue", "red"), alpha = 0.5)

# calculate z-score of vec_motif_data in distribution of vec_data elementwise
z_score <- (vec_motif_data - mean(vec_data)) / sd(vec_data)
mean(z_score)

z_score <- (vec_sub_data - mean(vec_data)) / sd(vec_data)
z_score
mean(z_score)
