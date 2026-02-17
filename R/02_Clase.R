# ==================================================
# 02-bioconductor-ejercicio.R - Código de la segunda clase
# Curso: RNA-Seq LCG-UNAM 2026
# Tema: Introducción a Bioconductor
# Exploración principal: ComplexHeatmap (por Gissel)
# ==================================================

## 1. INSTALACIÓN Y CONFIGURACIÓN BÁSICA ===================================

# Instalar BiocManager (solo una vez)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Verificar versión de Bioconductor
BiocManager::version()
# Debería mostrar: '3.22' (release actual)

# Instalar paquetes necesarios para la clase
paquetes_a_instalar <- c(
  "ComplexHeatmap", # Mi paquete: visualización avanzada
  "SummarizedExperiment", # Estructura de datos fundamental de BioC
  "circlize" # Para paletas de colores con ComplexHeatmap
)

# Instalar todos (ejecutar una vez)
BiocManager::install(paquetes_a_instalar)

## 2. EXPLORACIÓN DE BIOCONDUCTOR ==========================================

### 2.1 Información general de Bioconductor
# Ramas: release (3.22 estable) y devel (3.23 en desarrollo)
# Tipos de paquetes: Software, Annotation, Experiment Data, Workflows
# biocViews: sistema de clasificación jerárquica

### 2.2 SummarizedExperiment - Estructura de datos fundamental ------------
library(SummarizedExperiment)

# Crear datos de ejemplo: 200 genes, 6 muestras
n_genes <- 200
n_muestras <- 6
set.seed(42)

# Matriz de cuentas (simulando RNA-seq)
cuentas <- matrix(
  rnbinom(n_genes * n_muestras, size = 2, mu = 100),
  nrow = n_genes,
  ncol = n_muestras
)
rownames(cuentas) <- paste0("ENSG", sprintf("%08d", 1:n_genes))
colnames(cuentas) <- paste0("Sample", 1:n_muestras)

# Metadatos de las filas (genes)
row_data <- DataFrame(
  gene_id = rownames(cuentas),
  symbol = paste0("GENE", 1:n_genes),
  chromosome = sample(c("chr1", "chr2", "chr3"), n_genes, replace = TRUE)
)

# Metadatos de las columnas (muestras)
col_data <- DataFrame(
  sample_id = colnames(cuentas),
  condition = rep(c("Control", "Tratamiento"), each = 3),
  batch = rep(c(1, 2), 3),
  row.names = colnames(cuentas)
)

# Crear objeto SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = cuentas),
  rowData = row_data,
  colData = col_data
)

# Explorar el objeto
se
dim(se)
assayNames(se) # Nombres de los assays
colData(se) # Metadatos de muestras
rowData(se) # Metadatos de genes

# Acceder a los datos
head(assay(se, "counts")[, 1:4])

# Guardar objeto
saveRDS(se, file = "processed-data/SummarizedExperiment_ejemplo.rds")

## 3. EXPLORACIÓN PRINCIPAL: COMPLEXHEATMAP ================================
# ==================================================
# 02-bioconductor-ejercicio.R
# Tema: Introducción a Bioconductor con ComplexHeatmap
# ==================================================

# Verificar cita oficial del paquete principal
citation("ComplexHeatmap")

# También puedes verla en formato BibTeX si la necesitas para un paper
# print(citation("ComplexHeatmap"), bibtex = TRUE)

# Paquete explorado por: Gissel yoooooo
# Versión: 2.26.1 (BioC 3.22)
# Autor: Zuguang Gu
# Licencia: MIT + file LICENSE
# Propósito: Crear heatmaps complejos y altamente personalizables

library(ComplexHeatmap)
library(circlize) # Para paletas de colores

### 3.1 EJEMPLO BÁSICO: Heatmap simple ------------------------------------

# Crear matriz de ejemplo (15 genes, 8 muestras)
set.seed(123)
matriz_basica <- matrix(rnorm(120), nrow = 15, ncol = 8)
rownames(matriz_basica) <- paste0("Gene", 1:15)
colnames(matriz_basica) <- paste0("Muestra", 1:8)

# Heatmap básico
Heatmap(
  matriz_basica,
  name = "Expresión", # Título de la leyenda
  show_row_names = TRUE,
  show_column_names = TRUE
)

# Guardar
pdf("figuras/heatmap_basico.pdf", width = 8, height = 6)
Heatmap(matriz_basica, name = "Expresión")
dev.off()

print("Ejemplo básico completado: heatmap simple con 15 genes y 8 muestras")

### 3.2 EJEMPLO INTERMEDIO: Heatmap con anotaciones -----------------------

# Crear matriz más grande (30 genes, 12 muestras)
set.seed(456)
matriz_intermedia <- matrix(rnorm(360), nrow = 30, ncol = 12)
rownames(matriz_intermedia) <- paste0("Gene", 1:30)
colnames(matriz_intermedia) <- paste0("Muestra", 1:12)

# Crear anotaciones para las columnas (muestras)
anotacion_columnas <- HeatmapAnnotation(
  # Anotación categórica
  Condicion = c(rep("Control", 4), rep("Tratamiento", 4), rep("Vehiculo", 4)),
  # Anotación numérica
  Tiempo = rep(c(0, 24, 48), each = 4),
  # Definir colores
  col = list(
    Condicion = c(
      "Control" = "#1f78b4",
      "Tratamiento" = "#e31a23",
      "Vehiculo" = "#33a02c"
    ),
    Tiempo = colorRamp2(c(0, 24, 48), c("white", "yellow", "orange"))
  )
)

# Crear anotaciones para las filas (genes)
set.seed(789)
anotacion_filas <- rowAnnotation(
  # Anotaciones
  Regulacion = sample(c("Up", "Down", "Neutro"), 30, replace = TRUE),
  LogFC = runif(30, -2, 2),
  col = list(
    Regulacion = c("Up" = "red", "Down" = "blue", "Neutro" = "gray"),
    LogFC = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  )
)

# Heatmap con anotaciones
ht_intermedio <- Heatmap(
  matriz_intermedia,
  name = "Expresión",

  # Anotaciones
  top_annotation = anotacion_columnas,
  left_annotation = anotacion_filas,

  # Clustering
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",

  # Apariencia
  show_row_names = TRUE,
  show_column_names = TRUE,

  # Colores para el heatmap
  col = colorRamp2(c(-2, 0, 2), c("purple", "black", "yellow"))
)

# Dibujar
draw(ht_intermedio)

# Guardar
pdf("figuras/heatmap_intermedio_anotaciones.pdf", width = 12, height = 8)
draw(ht_intermedio)
dev.off()

print(
  "Ejemplo intermedio completado: heatmap con anotaciones de muestras y genes"
)

## 4. RESUMEN DE MI EXPLORACIÓN: COMPLEXHEATMAP ===========================

# ComplexHeatmap (v2.26.1) permite:
# 1. Heatmaps simples (ejemplo básico)
# 2. Anotaciones múltiples en columnas y filas (ejemplo intermedio)
# 3. Clustering personalizable
# 4. Control total de colores y estilos

# Estructura de objetos:
# - Heatmap-class: heatmap individual
# - HeatmapAnnotation-class: anotaciones para filas/columnas

print("=========================================")
print("FIN DEL SCRIPT - CLASE 2 COMPLETADA")
print("Exploración principal: ComplexHeatmap")
print("- Ejemplo básico: heatmap simple")
print("- Ejemplo intermedio: heatmap con anotaciones")
print("=========================================")

## 5. INFORMACIÓN DE LA SESIÓN ============================================

# Guardar información de la sesión
sink("processed-data/sesion_bioconductor_complexheatmap.txt")
print("===== VERSIÓN DE R =====")
version
print("===== PAQUETES CARGADOS =====")
sessionInfo()
print("===== VERSIÓN DE COMPLEXHEATMAP =====")
packageVersion("ComplexHeatmap")
sink()

# Guardar todo el entorno
save.image(file = "processed-data/clase2_complexheatmap.RData")

print("Todos los archivos guardados en:")
print("- figuras/heatmap_basico.pdf")
print("- figuras/heatmap_intermedio_anotaciones.pdf")
print("- processed-data/sesion_bioconductor_complexheatmap.txt")
print("- processed-data/clase2_complexheatmap.RData")
