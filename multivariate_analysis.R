source("scripts/setup.R")

# Cargar datos normalizados
norm_counts <- as.matrix(read.table("data/processed/normalized_counts.tsv", 
                                    header = TRUE, row.names = 1))
metadata <- read.csv("data/metadata.csv")
rownames(metadata) <- metadata$sample

# 1. PCA con todos los genes
pca <- prcomp(t(norm_counts), center = TRUE, scale. = TRUE)
pca_summary <- summary(pca)

# 2. Preparar datos para gráficos
pca_df <- as.data.frame(pca$x[, 1:5])
pca_df$sample <- rownames(pca_df)
pca_df <- merge(pca_df, metadata, by.x = "sample", by.y = "sample")

# 3. Gráfico PCA 2D
pdf("results/plots/PCA_analysis.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))

# PCA por estímulo
plot(pca_df$PC1, pca_df$PC2,
     main = paste0("PCA: PC1 (", round(pca_summary$importance[2,1]*100, 1), 
                   "%) vs PC2 (", round(pca_summary$importance[2,2]*100, 1), "%)"),
     xlab = paste0("PC1 (", round(pca_summary$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(pca_summary$importance[2,2]*100, 1), "%)"),
     pch = 19, cex = 1.5,
     col = stimulus_colors[pca_df$stimulus])
legend("topright", legend = names(stimulus_colors),
       col = stimulus_colors, pch = 19, title = "Estímulo")
text(pca_df$PC1, pca_df$PC2, labels = pca_df$time,
     pos = 3, cex = 0.7)

# PCA por tiempo
plot(pca_df$PC1, pca_df$PC2,
     main = paste0("PCA: PC1 (", round(pca_summary$importance[2,1]*100, 1), 
                   "%) vs PC2 (", round(pca_summary$importance[2,2]*100, 1), "%)"),
     xlab = paste0("PC1 (", round(pca_summary$importance[2,1]*100, 1), "%)"),
     ylab = paste0("PC2 (", round(pca_summary$importance[2,2]*100, 1), "%)"),
     pch = 19, cex = 1.5,
     col = time_colors[pca_df$time])
legend("topright", legend = names(time_colors),
       col = time_colors, pch = 19, title = "Tiempo")
text(pca_df$PC1, pca_df$PC2, labels = pca_df$stimulus,
     pos = 3, cex = 0.7)

# Varianza explicada
barplot(pca_summary$importance[2, 1:10] * 100,
        main = "Varianza explicada por componentes",
        xlab = "Componente principal",
        ylab = "Varianza explicada (%)",
        col = "steelblue", las = 2)

# Heatmap de las muestras (top 1000 genes variables)
var_genes <- apply(norm_counts, 1, var)
top_var_genes <- names(sort(var_genes, decreasing = TRUE))[1:1000]
top_var_matrix <- norm_counts[top_var_genes, ]

annotation_df <- data.frame(
  Stimulus = metadata$stimulus,
  Time = metadata$time,
  row.names = metadata$sample
)

pheatmap(top_var_matrix,
         main = "Heatmap: Top 1000 genes más variables",
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         show_rownames = FALSE,
         show_colnames = TRUE,
         annotation_col = annotation_df,
         annotation_colors = list(
           Stimulus = stimulus_colors,
           Time = time_colors
         ),
         color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(100))

dev.off()

# 4. Guardar resultados del PCA
write.csv(pca_summary$importance, "results/tables/PCA_variance.csv")
write.csv(pca_df, "results/tables/PCA_coordinates.csv", row.names = FALSE)

print("Análisis multivariante completado")
