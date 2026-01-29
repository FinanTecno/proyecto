source("scripts/setup.R")

# Cargar resultados DE
de_results <- read.csv("results/tables/Differential_expression_all_comparisons.csv")

# Función para análisis de enriquecimiento
run_enrichment <- function(gene_list, comparison_name, direction) {
  cat("\nEnriquecimiento para:", comparison_name, "-", direction, "\n")
  
  # Convertir IDs (en datos reales, necesitarías mapear ENSG a Entrez)
  # Para este ejemplo, usaré IDs simulados
  
  if (length(gene_list) < 10) {
    cat("Muy pocos genes para enriquecimiento\n")
    return(NULL)
  }
  
  # GO Biological Process
  ego_bp <- enrichGO(gene = gene_list,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENSEMBL",  # Cambiar según tus IDs
                     ont = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.1,
                     readable = TRUE)
  
  # KEGG pathways
  ekegg <- enrichKEGG(gene = gene_list,
                      organism = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05)
  
  # Guardar resultados
  if (!is.null(ego_bp) && nrow(ego_bp) > 0) {
    write.csv(as.data.frame(ego_bp),
              paste0("results/enrichment/", 
                     gsub(" ", "_", comparison_name), "_",
                     direction, "_GO_BP.csv"),
              row.names = FALSE)
  }
  
  if (!is.null(ekegg) && nrow(ekegg) > 0) {
    write.csv(as.data.frame(ekegg),
              paste0("results/enrichment/",
                     gsub(" ", "_", comparison_name), "_",
                     direction, "_KEGG.csv"),
              row.names = FALSE)
  }
  
  return(list(GO_BP = ego_bp, KEGG = ekegg))
}

# Realizar enriquecimiento para cada comparación significativa
for (comp in unique(de_results$comparison)) {
  comp_data <- de_results %>% filter(comparison == comp)
  
  # Genes up-regulados
  up_genes <- comp_data %>% 
    filter(direction == "Up") %>% 
    pull(gene)
  
  # Genes down-regulados
  down_genes <- comp_data %>% 
    filter(direction == "Down") %>% 
    pull(gene)
  
  if (length(up_genes) > 0) {
    run_enrichment(up_genes, comp, "Up")
  }
  
  if (length(down_genes) > 0) {
    run_enrichment(down_genes, comp, "Down")
  }
}

# 2. Heatmap de los top DEGs
norm_counts <- as.matrix(read.table("data/processed/normalized_counts.tsv", 
                                    header = TRUE, row.names = 1))

# Seleccionar top 50 genes por comparación
top_genes <- de_results %>%
  filter(significant) %>%
  group_by(comparison) %>%
  arrange(padj) %>%
  slice_head(n = 12) %>%
  pull(gene) %>%
  unique()

if (length(top_genes) > 0) {
  heatmap_data <- norm_counts[top_genes, ]
  
  # Anotaciones
  annotation_df <- data.frame(
    Stimulus = metadata$stimulus,
    Time = metadata$time,
    row.names = metadata$sample
  )
  
  pdf("results/plots/Top_DEGs_heatmap.pdf", width = 12, height = 10)
  pheatmap(heatmap_data,
           main = "Top genes diferencialmente expresados",
           scale = "row",
           clustering_distance_rows = "euclidean",
           clustering_distance_cols = "euclidean",
           show_rownames = TRUE,
           show_colnames = TRUE,
           fontsize_row = 8,
           annotation_col = annotation_df,
           annotation_colors = list(
             Stimulus = stimulus_colors,
             Time = time_colors
           ),
           color = colorRampPalette(rev(brewer.pal(10, "RdBu")))(100))
  dev.off()
}

# 3. Gráfico de barras de número de DEGs
de_summary <- de_results %>%
  filter(direction %in% c("Up", "Down")) %>%
  group_by(comparison, stimulus, time_comparison, direction) %>%
  summarise(count = n(), .groups = "drop")

pdf("results/plots/DEGs_summary_barplot.pdf", width = 10, height = 6)
ggplot(de_summary, aes(x = comparison, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = count), 
            position = position_dodge(width = 0.9),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("Down" = "#0072B2", "Up" = "#D55E00")) +
  labs(title = "Número de genes diferencialmente expresados",
       x = "Comparación",
       y = "Número de genes") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

print("Análisis de enriquecimiento completado")
