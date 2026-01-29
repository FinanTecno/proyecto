source("scripts/setup.R")

# Cargar datos
expression_matrix <- as.matrix(read.table("data/raw/expression_matrix.tsv", 
                                          header = TRUE, row.names = 1))
metadata <- read.csv("data/metadata.csv")
rownames(metadata) <- metadata$sample
metadata <- metadata[colnames(expression_matrix), ]

# Definir las comparaciones
comparisons <- list(
  LPS_4h_vs_0h = list(stimulus = "LPS", time1 = "4h", time2 = "0h"),
  LPS_24h_vs_0h = list(stimulus = "LPS", time1 = "24h", time2 = "0h"),
  Bglucan_4h_vs_0h = list(stimulus = "Bglucan", time1 = "4h", time2 = "0h"),
  Bglucan_24h_vs_0h = list(stimulus = "Bglucan", time1 = "24h", time2 = "0h")
)

# Función para realizar DE
run_de_analysis <- function(comp_name, comp_info) {
  cat("\n=== Análisis:", comp_name, "===\n")
  
  # Filtrar muestras para esta comparación
  samples <- metadata %>%
    filter(stimulus == comp_info$stimulus,
           time %in% c(comp_info$time1, comp_info$time2))
  
  counts_subset <- expression_matrix[, samples$sample]
  
  # Crear objeto DESeq2
  colData <- samples %>%
    mutate(time_factor = factor(time, levels = c(comp_info$time2, comp_info$time1)))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_subset,
    colData = colData,
    design = ~ time_factor
  )
  
  # Filtrar y normalizar
  dds <- dds[rowSums(counts(dds) >= 10) >= 3, ]
  dds <- DESeq(dds)
  
  # Resultados
  res <- results(dds, alpha = 0.05, 
                 contrast = c("time_factor", comp_info$time1, comp_info$time2))
  
  # Ordenar por valor-p ajustado
  res <- res[order(res$padj), ]
  
  # Añadir significancia
  res_df <- as.data.frame(res) %>%
    mutate(
      gene = rownames(res),
      comparison = comp_name,
      stimulus = comp_info$stimulus,
      time_comparison = paste(comp_info$time1, "vs", comp_info$time2),
      significant = padj < 0.05 & abs(log2FoldChange) > 1,
      direction = case_when(
        significant & log2FoldChange > 0 ~ "Up",
        significant & log2FoldChange < 0 ~ "Down",
        TRUE ~ "Non-sig"
      )
    )
  
  return(res_df)
}

# Ejecutar todas las comparaciones
all_results <- list()
for (comp_name in names(comparisons)) {
  all_results[[comp_name]] <- run_de_analysis(comp_name, comparisons[[comp_name]])
}

# Combinar todos los resultados
combined_results <- do.call(rbind, all_results)
rownames(combined_results) <- NULL

# Guardar resultados
write.csv(combined_results, 
          "results/tables/Differential_expression_all_comparisons.csv",
          row.names = FALSE)

# Resumen por comparación
summary_de <- combined_results %>%
  group_by(comparison, stimulus, time_comparison, direction) %>%
  filter(direction %in% c("Up", "Down")) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = direction, values_from = count, values_fill = 0)

write.csv(summary_de, 
          "results/tables/DE_summary_counts.csv",
          row.names = FALSE)

# 7. Volcano plots para cada comparación
pdf("results/plots/Volcano_plots.pdf", width = 12, height = 10)
par(mfrow = c(2, 2))

for (comp in names(comparisons)) {
  comp_data <- all_results[[comp]]
  
  plot(comp_data$log2FoldChange, -log10(comp_data$pvalue),
       main = paste("Volcano Plot:", comp),
       xlab = "log2(Fold Change)",
       ylab = "-log10(p-value)",
       pch = 19, cex = 0.6,
       col = case_when(
         comp_data$direction == "Up" ~ "#D55E00",
         comp_data$direction == "Down" ~ "#0072B2",
         TRUE ~ "#CCCCCC"
       ))
  
  abline(v = c(-1, 1), lty = 2, col = "gray")
  abline(h = -log10(0.05), lty = 2, col = "gray")
  
  # Etiquetar top 10 genes
  top_genes <- comp_data %>%
    filter(significant) %>%
    arrange(padj) %>%
    head(10)
  
  if (nrow(top_genes) > 0) {
    text(top_genes$log2FoldChange, 
         -log10(top_genes$pvalue),
         labels = top_genes$gene,
         pos = 3, cex = 0.5)
  }
  
  legend("topright", 
         legend = c(paste("Up:", sum(comp_data$direction == "Up")),
                    paste("Down:", sum(comp_data$direction == "Down"))),
         col = c("#D55E00", "#0072B2"), pch = 19)
}

dev.off()

print("Análisis de expresión diferencial completado")
