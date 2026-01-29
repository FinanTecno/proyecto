source("scripts/setup.R")

# Método 1: Descargar desde GEO (si tienes acceso a internet)
# library(GEOquery)
# gse <- getGEO("GSE111003", destdir = "data/raw")
# expression_matrix <- exprs(gse[[1]])

# Método 2: Cargar matriz local (si ya la tienes)
# NOTA: Necesitas descargar manualmente la matriz de GSE111003
# y guardarla como data/raw/expression_matrix.tsv

# Para este ejemplo, crearé una matriz de ejemplo
# En la práctica real, usarías los datos reales

set.seed(123)
n_genes <- 20000
n_samples <- 18

# Crear nombres de genes (ejemplo)
gene_ids <- paste0("ENSG", formatC(1:n_genes, width = 11, flag = "0"))

# Crear matriz de expresión (simulación)
base_expression <- rnbinom(n_genes * n_samples, mu = 100, size = 10)
expression_matrix <- matrix(base_expression, nrow = n_genes, ncol = n_samples)

# Añadir efectos diferenciales
# Genes diferenciales para LPS 4h
de_genes_lps_4h <- sample(1:n_genes, 500)
expression_matrix[de_genes_lps_4h, 4:6] <- expression_matrix[de_genes_lps_4h, 4:6] * runif(500, 1.5, 4)

# Genes diferenciales para LPS 24h
de_genes_lps_24h <- sample(setdiff(1:n_genes, de_genes_lps_4h), 800)
expression_matrix[de_genes_lps_24h, 7:9] <- expression_matrix[de_genes_lps_24h, 7:9] * runif(800, 2, 6)

# Genes diferenciales para B-glucano
de_genes_bglucan <- sample(1:n_genes, 300)
expression_matrix[de_genes_bglucan, 13:15] <- expression_matrix[de_genes_bglucan, 13:15] * runif(300, 1.3, 3)

# Añadir nombres
rownames(expression_matrix) <- gene_ids
colnames(expression_matrix) <- c(
  "GSM3031851", "GSM3031852", "GSM3031853",
  "GSM3031854", "GSM3031855", "GSM3031856",
  "GSM3031857", "GSM3031858", "GSM3031859",
  "GSM3031860", "GSM3031861", "GSM3031862",
  "GSM3031863", "GSM3031864", "GSM3031865",
  "GSM3031866", "GSM3031867", "GSM3031868"
)

# Guardar matriz
write.table(expression_matrix, 
            "data/raw/expression_matrix.tsv",
            sep = "\t", quote = FALSE)

# Cargar metadatos
metadata <- read.csv("data/metadata.csv", stringsAsFactors = FALSE)
metadata$condition <- factor(metadata$condition, 
                             levels = unique(metadata$condition))

print("Datos cargados exitosamente")
print(paste("Dimensiones de la matriz:", 
            nrow(expression_matrix), "genes x", 
            ncol(expression_matrix), "muestras"))
