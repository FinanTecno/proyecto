source("scripts/setup.R")

# Cargar datos
expression_matrix <- as.matrix(read.table("data/raw/expression_matrix.tsv", 
                                          header = TRUE, row.names = 1))
metadata <- read.csv("data/metadata.csv")
rownames(metadata) <- metadata$sample

# Asegurar que el orden coincide
metadata <- metadata[colnames(expression_matrix), ]

# 1. Crear objeto DESeq2
dds <- DESeqDataSetFromMatrix(
  countData = expression_matrix,
  colData = metadata,
  design = ~ stimulus + time_numeric + stimulus:time_numeric
)

# 2. Filtrar genes de baja expresión
# Mantener genes con al menos 10 cuentas en al menos 3 muestras
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
print(paste("Genes después del filtrado:", nrow(dds)))

# 3. Normalización con DESeq2
dds <- DESeq(dds)

# 4. Obtener datos normalizados
vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
norm_counts <- assay(vsd)

# 5. Guardar datos normalizados
write.table(norm_counts, 
            "data/processed/normalized_counts.tsv",
            sep = "\t", quote = FALSE)

# 6. Boxplot de datos normalizados
pdf("results/plots/normalization_boxplot.pdf", width = 12, height = 6)
boxplot(norm_counts, 
        main = "Distribución después de normalización (VST)",
        xlab = "Muestras", ylab = "Valores transformados",
        las = 2, col = stimulus_colors[metadata$stimulus],
        cex.axis = 0.7)
legend("topright", legend = names(stimulus_colors),
       fill = stimulus_colors, title = "Estímulo")
dev.off()

# 7. Gráfico de dispersión de medias vs varianzas
pdf("results/plots/mean_variance_plot.pdf", width = 8, height = 6)
mean_counts <- rowMeans(counts(dds, normalized = TRUE))
var_counts <- apply(counts(dds, normalized = TRUE), 1, var)
plot(log2(mean_counts + 1), log2(var_counts + 1),
     main = "Media vs Varianza (log2)",
     xlab = "Media de cuentas (log2)",
     ylab = "Varianza (log2)",
     pch = 16, cex = 0.5, col = rgb(0, 0, 1, 0.3))
abline(a = 0, b = 1, col = "red", lwd = 2)
legend("topleft", legend = "y = x (Poisson)", col = "red", lwd = 2)
dev.off()

print("Normalización completada")
