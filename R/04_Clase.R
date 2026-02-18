# ==================================================
# 04-recount3.R - Código de la cuarta clase
# Curso: RNA-Seq LCG-UNAM 2026
# Tema: Datos de RNA-seq a través de recount3
# ==================================================

## 1. INSTALACIÓN Y CARGA ==================================================

# Instalar recount3 (solo una vez)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("recount3")

# Cargar paquetes necesarios
library("recount3")
library("iSEE")
library("ggplot2")
library("sessioninfo")

## 2. EXPLORAR PROYECTOS DISPONIBLES =======================================

# Revisar todos los proyectos con datos de humano
human_projects <- available_projects()

# Ver estructura
dim(human_projects)
head(human_projects)

# Opcional: proyectos de ratón
# mouse_projects <- available_projects(organism = "mouse")

## 3. SELECCIONAR UN PROYECTO ESPECÍFICO ===================================
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("interactiveDisplayBase")
# Encontrar el proyecto SRP009615
proj_info <- subset(
  human_projects,
  project == "SRP009615" & project_type == "data_sources"
)
proj_info

# Opción interactiva (seleccionar manualmente)
if (interactive()) {
  proj_info_interactive <- interactiveDisplayBase::display(human_projects)
  # Verificar que solo seleccionaste un renglón
  stopifnot(nrow(proj_info_interactive) == 1)
}

## 4. DESCARGAR DATOS CON CREATE_RSE() =====================================

# Crear objeto RangedSummarizedExperiment
rse_gene_SRP009615 <- create_rse(proj_info)

# Explorar el objeto
rse_gene_SRP009615
dim(rse_gene_SRP009615)
assayNames(rse_gene_SRP009615) # "raw_counts"

## 5. TRANSFORMAR CUENTAS ==================================================

# Convertir cuentas por nucleótido a cuentas por lectura
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)

# Verificar que ahora tenemos dos assays
assayNames(rse_gene_SRP009615)
# [1] "raw_counts" "counts"

# Ver primeras filas de las cuentas transformadas
head(assay(rse_gene_SRP009615, "counts"))

## 6. EXPANDIR METADATOS DE SRA ============================================

# Expandir atributos para hacerlos más legibles
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)

# Ver los atributos expandidos
colData(rse_gene_SRP009615)[,
  grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

## 7. EXPLORAR EL OBJETO COMPLETO ==========================================

# Información de los genes
rowRanges(rse_gene_SRP009615)

# Metadatos de las muestras
colData(rse_gene_SRP009615)

# Nombres completos de las columnas en colData
colnames(colData(rse_gene_SRP009615))

## 8. EJERCICIO: REPRODUCIR IMAGEN CON iSEE ================================

### 8.1 Preparar datos para iSEE (simplificar nombres) --------------------

# Crear versiones simplificadas de las variables
colData(rse_gene_SRP009615)$treatment <- colData(
  rse_gene_SRP009615
)$sra_attribute.treatment
colData(rse_gene_SRP009615)$shRNA <- colData(
  rse_gene_SRP009615
)$sra_attribute.shRNA_expression
colData(rse_gene_SRP009615)$cells <- colData(
  rse_gene_SRP009615
)$sra_attribute.cells
colData(rse_gene_SRP009615)$source <- colData(
  rse_gene_SRP009615
)$sra_attribute.source_name

# Verificar nuevas variables
head(colData(rse_gene_SRP009615)[, c("treatment", "shRNA", "cells", "source")])

### 8.2 Lanzar iSEE para exploración interactiva --------------------------

if (interactive()) {
  iSEE::iSEE(rse_gene_SRP009615)
}

### 8.3 Pasos dentro de iSEE para reproducir la imagen --------------------
#
# 1. Seleccionar "Column data plot" como tipo de panel
# 2. Configurar eje X: "Column data" -> "treatment"
# 3. Configurar eje Y: "Feature assay" -> seleccionar un gen (ej. el primero)
# 4. Configurar color: "Column data" -> "shRNA"
# 5. Ajustar tipo de gráfico: "Boxplot" o "Violin"
# 6. Guardar como PDF con el botón "Download"

### 8.4 Explorar otros genes ----------------------------------------------

# Función para graficar cualquier gen
plot_gene_expression <- function(rse, gene_id) {
  if (!gene_id %in% rownames(rse)) {
    stop("Gen no encontrado")
  }

  exp_data <- data.frame(
    expression = as.numeric(assay(rse, "counts")[gene_id, ]),
    treatment = colData(rse)$treatment,
    shRNA = colData(rse)$shRNA
  )

  ggplot(exp_data, aes(x = treatment, y = expression, fill = shRNA)) +
    geom_boxplot(alpha = 0.7) +
    geom_jitter(position = position_jitterdodge(), size = 1.5, alpha = 0.5) +
    theme_bw() +
    labs(
      title = paste("Expresión de", gene_id),
      x = "Tratamiento",
      y = "Cuentas por lectura",
      fill = "shRNA"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Probar con algunos genes (primeros 5 genes)
for (i in 1:5) {
  gen <- rownames(rse_gene_SRP009615)[i]
  print(plot_gene_expression(rse_gene_SRP009615, gen))
  ggsave(paste0("figuras/boxplot_gen_", i, ".pdf"), width = 8, height = 5)
}

## 9. ANÁLISIS ADICIONAL: RESUMEN POR CONDICIÓN ===========================

# Resumen estadístico por tratamiento y shRNA
library("dplyr")

# Crear data frame con todos los genes (primeros 100 para ejemplo)
genes_sample <- rownames(rse_gene_SRP009615)[1:100]
expresion_matrix <- t(assay(rse_gene_SRP009615, "counts")[genes_sample, ])

df_summary <- as.data.frame(expresion_matrix)
df_summary$treatment <- colData(rse_gene_SRP009615)$treatment
df_summary$shRNA <- colData(rse_gene_SRP009615)$shRNA

# Media de expresión por grupo para el primer gen
df_summary %>%
  group_by(treatment, shRNA) %>%
  summarise(
    mean_expr = mean(V1), # V1 es el primer gen
    sd_expr = sd(V1),
    n = n()
  )

## 10. GUARDAR OBJETOS Y SESIÓN ============================================

# Guardar el objeto RSE para uso futuro
dir.create("processed-data", showWarnings = FALSE)
saveRDS(rse_gene_SRP009615, file = "processed-data/rse_SRP009615.rds")

# Guardar información de la sesión
sink("processed-data/sesion_clase4.txt")
print("===== VERSIÓN DE R =====")
version
print("===== PAQUETES CARGADOS =====")
sessionInfo()
print("===== PROYECTO DESCARAGDO =====")
proj_info
print("===== DIMENSIONES DEL OBJETO =====")
dim(rse_gene_SRP009615)
sink()

# Guardar imagen completa del entorno
save.image(file = "processed-data/clase4_recount3.RData")

## 11. VERIFICACIONES FINALES ==============================================

# Verificar que los archivos se guardaron correctamente
cat("\n=========================================\n")
cat("ARCHIVOS GENERADOS:\n")
cat("=========================================\n")
cat("- processed-data/rse_SRP009615.rds\n")
cat("- processed-data/sesion_clase4.txt\n")
cat("- processed-data/clase4_recount3.RData\n")
cat("- figuras/boxplot_recount3_ejercicio.pdf\n")
cat("- figuras/boxplot_recount3_ejercicio.png\n")
for (i in 1:5) {
  cat(paste0("- figuras/boxplot_gen_", i, ".pdf\n"))
}
cat("=========================================\n")
cat("FIN DEL SCRIPT - CLASE 4 COMPLETADA\n")
cat("=========================================\n")
