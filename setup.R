# Instalar paquetes necesarios (ejecutar solo una vez)
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install(c(
#   "DESeq2", "edgeR", "limma",
#   "ggplot2", "pheatmap", "RColorBrewer",
#   "dplyr", "tidyr", "reshape2",
#   "clusterProfiler", "org.Hs.eg.db",
#   "GEOquery", "ggrepel", "factoextra"
# ))

# Cargar librerías
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(factoextra)
library(reshape2)

# Configurar tema para gráficos
theme_set(theme_minimal(base_size = 12) +
            theme(plot.title = element_text(hjust = 0.5, face = "bold"),
                  legend.position = "bottom"))

# Crear paleta de colores
stimulus_colors <- c("LPS" = "#E41A1C", "Bglucan" = "#377EB8")
time_colors <- c("0h" = "#F0E442", "4h" = "#E69F00", "24h" = "#D55E00")

# Directorios
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/enrichment", recursive = TRUE, showWarnings = FALSE)
