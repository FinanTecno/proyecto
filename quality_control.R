source("scripts/setup.R")

# Cargar datos
expression_matrix <- as.matrix(read.table("data/raw/expression_matrix.tsv", 
                                          header = TRUE, row.names = 1))
metadata <- read.csv("data/metadata.csv")

# 1. Distribución de cuentas por muestra (raw)
counts_per_sample <- colSums(expression_matrix)
genes_per_sample <- colSums(expression_matrix > 0)

qc_summary <- data.frame(
  Sample = colnames(expression_matrix),
  Total_Counts = counts_per_sample,
  Genes_Detected = genes_per_sample
)

qc_summary <- merge(qc_summary, metadata, by.x = "Sample", by.y = "sample")

# Guardar resumen
write.csv(qc_summary, "results/tables/QC_summary.csv", row.names = FALSE)

# 2. Boxplot de cuentas raw
pdf("results/plots/QC_boxplot_raw.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))

# Boxplot log10
boxplot(log10(expression_matrix + 1), 
        main = "Distribución de cuentas (log10) - RAW",
        xlab = "Muestras", ylab = "log10(counts + 1)",
        las = 2, col = stimulus_colors[qc_summary$stimulus],
        cex.axis = 0.7)
legend("topright", legend = names(stimulus_colors),
       fill = stimulus_colors, title = "Estímulo")

# Genes detectados
barplot(qc_summary$Genes_Detected/1000, 
        names.arg = qc_summary$Sample,
        main = "Genes detectados por muestra",
        xlab = "Muestras", ylab = "Número de genes (x1000)",
        las = 2, col = time_colors[qc_summary$time],
        cex.names = 0.7)
legend("topright", legend = names(time_colors),
       fill = time_colors, title = "Tiempo")
dev.off()

# 3. Matriz de correlación
cor_matrix <- cor(expression_matrix, method = "spearman")

pdf("results/plots/QC_correlation_heatmap.pdf", width = 10, height = 8)
pheatmap(cor_matrix,
         main = "Correlación entre muestras (Spearman)",
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = FALSE,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_row = data.frame(
           Stimulus = metadata$stimulus,
           Time = metadata$time,
           row.names = metadata$sample
         ),
         annotation_colors = list(
           Stimulus = stimulus_colors,
           Time = time_colors
         ))
dev.off()

print("Control de calidad completado")
