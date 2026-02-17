# ==================================================
# 03-objetos-bioconductor.R - Código de la tercera clase
# Curso: RNA-Seq LCG-UNAM 2026
# Tema: Objetos de Bioconductor (SummarizedExperiment, iSEE, spatialLIBD)
# ==================================================

## 1. SUMMARIZEDEXPERIMENT: CREACIÓN Y EXPLORACIÓN =========================

# Cargar librerías necesarias
library("SummarizedExperiment")
library("GenomicRanges")

### 1.1 Crear un objeto SummarizedExperiment desde cero --------------------

# Definir dimensiones
nrows <- 200 # genes
ncols <- 6 # muestras
set.seed(20210223)

# Matriz de cuentas (assay)
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rownames(counts) <- paste0("gene_", 1:nrows)
colnames(counts) <- LETTERS[1:ncols]

# Información genómica de los genes (rowRanges)
rowRanges <- GRanges(
  seqnames = rep(c("chr1", "chr2"), c(50, 150)),
  ranges = IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
  strand = sample(c("+", "-"), 200, TRUE),
  feature_id = sprintf("ID%03d", 1:200)
)
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

# Metadatos de las muestras (colData)
colData <- DataFrame(
  Treatment = rep(c("ChIP", "Input"), 3),
  row.names = LETTERS[1:6]
)

# Ensamblar todo en un objeto SummarizedExperiment
rse <- SummarizedExperiment(
  assays = SimpleList(counts = counts),
  rowRanges = rowRanges,
  colData = colData
)

### 1.2 Explorar el objeto -----------------------------------------------

# Vista general
rse

# Dimensiones
dim(rse)

# Nombres de filas y columnas
dimnames(rse)

# Nombres de los assays disponibles
assayNames(rse)

# Ver primeras filas de la matriz de cuentas
head(assay(rse, "counts"))

# Información genómica completa
rowRanges(rse)

# Metadatos de los genes (tabla plana)
rowData(rse)

# Metadatos de las muestras
colData(rse)

## 2. EJERCICIO: MANIPULACIÓN DE SUMMARIZEDEXPERIMENT =====================

# Comando 1: Seleccionar primeros 2 genes, todas las muestras
rse[1:2, ]
# Resultado: Objeto con genes gene_1 y gene_2 para las 6 muestras

# Comando 2: Seleccionar todos los genes, solo muestras A, D y F
rse[, c("A", "D", "F")]
# Resultado: Objeto con 200 genes para las muestras A, D y F

## 3. EXPLORACIÓN INTERACTIVA CON iSEE ====================================

# Cargar iSEE
library("iSEE")

# Lanzar app interactiva con nuestro objeto rse
# (Solo funciona en sesión interactiva)
if (interactive()) {
  iSEE::iSEE(rse)
}

## 4. EJERCICIO PRÁCTICO CON spatialLIBD ===================================

# Cargar librerías adicionales
library("spatialLIBD")
library("lobstr") # Para ver tamaño de objetos
library("scater") # Para plotting de SingleCellExperiment

### 4.1 Descargar y explorar datos de spatialLIBD -------------------------

# Descargar datos (se guardan en caché automáticamente)
sce_layer <- spatialLIBD::fetch_data("sce_layer")
# Mensaje: adding rname 'https://www.dropbox.com/...'

# Explorar objeto descargado
sce_layer
# class: SingleCellExperiment
# dim: 22331 76
# assays(2): counts logcounts
# rowData names(10): source type ...
# colData names(13): sample_name layer_guess ...

# Ver tamaño en memoria
lobstr::obj_size(sce_layer)
# ~33.99 MB

### 4.2 Lanzar iSEE para exploración interactiva -------------------------

if (interactive()) {
  iSEE::iSEE(sce_layer)
}

### 4.3 Resolver preguntas del ejercicio (solo código) ----------------

# Asegurar que tenemos logcounts
if (!"logcounts" %in% assayNames(sce_layer)) {
  sce_layer <- scater::logNormCounts(sce_layer)
}

# IDs de los genes de interés
id_MOBP <- "ENSG00000168314" # MOBP
id_MBP <- "ENSG00000197971" # MBP
id_PCP4 <- "ENSG00000183036" # PCP4

#### Pregunta 1: Heatmap de MOBP, MBP y PCP4 -----------------------------
# Extraer matriz de expresión para los 3 genes
genes_interes <- c(id_MOBP, id_MBP, id_PCP4)
expr_matrix <- assay(sce_layer[genes_interes, ], "logcounts")
rownames(expr_matrix) <- c("MOBP", "MBP", "PCP4") # Renombrar para claridad

# Crear heatmap básico
heatmap(
  expr_matrix,
  main = "Patrón de expresión: MOBP, MBP, PCP4",
  xlab = "Muestras",
  ylab = "Genes"
)

# Conclusión: MOBP y MBP tienen patrones similares (se agrupan juntos)
# PCP4 tiene patrón diferente

#### Pregunta 2: ¿En qué capas se expresan más MOBP y MBP? --------------

# Boxplot para MOBP por capa
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("scater")

library(scater)


plotExpression(
  sce_layer,
  features = id_MOBP,
  x = "layer_guess_reordered_short",
  colour_by = "layer_guess_reordered_short"
) +
  ggtitle("Expresión de MOBP por capa cortical")

# Boxplot para MBP por capa
plotExpression(
  sce_layer,
  features = id_MBP,
  x = "layer_guess_reordered_short",
  colour_by = "layer_guess_reordered_short"
) +
  ggtitle("Expresión de MBP por capa cortical")

# Conclusión: MOBP y MBP se expresan más en:
# - Materia blanca (WM)
# - Capa 6 (Layer6)

### 4.4 Verificar con estadísticas resumidas -----------------------------

# Obtener información de capas
capas <- colData(sce_layer)$layer_guess_reordered_short

# Calcular media de expresión por capa para MOBP
tapply(assay(sce_layer[id_MOBP, ], "logcounts")[1, ], capas, mean)

# Calcular media de expresión por capa para MBP
tapply(assay(sce_layer[id_MBP, ], "logcounts")[1, ], capas, mean)

## 5. INFORMACIÓN DE LA SESIÓN ============================================

# Guardar información de la sesión para reproducibilidad
sink("processed-data/sesion_clase3.txt")
print("===== VERSIÓN DE R =====")
version
print("===== PAQUETES CARGADOS =====")
sessionInfo()
sink()

# Guardar objetos creados
saveRDS(rse, file = "processed-data/rse_creado_en_clase.rds")
saveRDS(sce_layer, file = "processed-data/sce_layer_spatialLIBD.rds")

# Guardar imagen completa del entorno
save.image(file = "processed-data/clase3_completa.RData")

print("=========================================")
print("FIN DEL SCRIPT - CLASE 3 COMPLETADA")
print("Objetos creados y ejercicios resueltos")
print("=========================================")
