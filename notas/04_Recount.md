# Notas de Clase: Datos de RNA-seq a través de recount3

**Curso:** RNA-Seq LCG-UNAM 2026

## Objetivos de la Clase de Hoy

1.  Comprender la filosofía de `recount3`: democratizar el acceso a datos públicos de RNA-seq.
2.  Aprender a buscar, identificar y descargar proyectos de interés usando `recount3`.
3.  Explorar la estructura de los objetos `RangedSummarizedExperiment` generados por `recount3`.
4.  Realizar transformaciones básicas de los datos (de cuentas por nucleótido a cuentas por lectura).
5.  Utilizar `iSEE` para visualizar y explorar los datos descargados, reproduciendo una figura específica.

---

## 1. Introducción a recount3

`recount3` es un recurso masivo que proporciona datos de RNA-seq **uniformemente procesados** para más de **700,000 muestras** de humano y ratón. Su objetivo es democratizar el acceso a los datos, permitiendo que cualquier persona con R pueda analizarlos sin necesidad de infraestructura de cómputo de alto rendimiento (HPC).

### 1.1 Antecedentes: La familia recount

| Proyecto | Muestras | Año | Características |
|:---|:---:|:---:|:---|
| **ReCount** | ~20 estudios | 2011 | Primer esfuerzo, datos limitados. |
| **recount** | ~70,000 | 2017 | Procesamiento uniforme de estudios públicos. |
| **recount3** | **>700,000** | 2021 | Humano y ratón, múltiples fuentes (SRA, GTEx, TCGA). |

### 1.2 Publicaciones clave

*   **Artículo principal recount3:** Wilks C, et al. (2021). "recount3: summaries and queries for large-scale RNA-seq expression and splicing". *Genome Biology*. [DOI: 10.1186/s13059-021-02533-6](https://doi.org/10.1186/s13059-021-02533-6)
*   **Pre-print:** mayo 2021. [DOI: 10.1101/2021.05.21.445138](https://doi.org/10.1101/2021.05.21.445138)
*   **Artículo recount (2017):** Collado-Torres L, et al. (2017). "Reproducible RNA-seq analysis using recount2". *Nature Biotechnology*. [DOI: 10.1038/nbt.3838](https://doi.org/10.1038/nbt.3838)

---

## 2. Uso de recount3 en R

### 2.1 Instalación y carga

```r
# Instalar recount3 (si no está instalado)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("recount3")

# Cargar el paquete
library("recount3")
```

### 2.2 Identificar proyectos disponibles

El primer paso es explorar qué proyectos están disponibles. Usamos la función `available_projects()`.

```r
# Revisar todos los proyectos con datos de humano
human_projects <- available_projects()

# Ver estructura del data frame
dim(human_projects)
head(human_projects)

# Podemos explorar proyectos de ratón también
# mouse_projects <- available_projects(organism = "mouse")
```

### 2.3 Seleccionar un proyecto de interés

Podemos seleccionar un proyecto específico, por ejemplo, el estudio **SRP009615** (un estudio de RNA-seq en células K562 con shRNA).

```r
# Encontrar el proyecto específico
proj_info <- subset(
    human_projects,
    project == "SRP009615" & project_type == "data_sources"
)
proj_info
```

**Selección interactiva (opcional):** También podemos usar una interfaz interactiva para elegir el proyecto.

```r
# Explorar proyectos de forma interactiva
if(interactive()) {
    proj_info_interactive <- interactiveDisplayBase::display(human_projects)
    # Seleccionar un renglón y hacer click en "send"
    stopifnot(nrow(proj_info_interactive) == 1)  # Verificar que solo uno está seleccionado
}
```

### 2.4 Descargar los datos con `create_rse()`

Una vez identificado el proyecto, descargamos los datos usando `create_rse()`. Esto crea un objeto `RangedSummarizedExperiment` (RSE) con los datos a nivel de genes.

```r
# Crear objeto RSE con datos a nivel de genes
rse_gene_SRP009615 <- create_rse(proj_info)

# Explorar el objeto
rse_gene_SRP009615
# class: RangedSummarizedExperiment 
# dim: 63856 12 
# assays(1): raw_counts
# rownames(63856): ENSG00000278704.1 ENSG00000277400.1 ...
# colnames(12): SRR387777 SRR387778 ... SRR389077 SRR389078
# colData names(175): rail_id external_id ... BigWigURL
```

### 2.5 Transformar las cuentas

Los datos en `recount3` vienen como **cuentas por nucleótido** (suma de lecturas por base). Para la mayoría de los análisis, necesitamos **cuentas por lectura** (similar a lo que produce featureCounts o HTSeq). Usamos `compute_read_counts()` para esta transformación.

```r
# Convertir cuentas por nucleótido a cuentas por lectura
assay(rse_gene_SRP009615, "counts") <- compute_read_counts(rse_gene_SRP009615)

# Verificar que ahora tenemos dos assays
assayNames(rse_gene_SRP009615)
# [1] "raw_counts" "counts"
```

### 2.6 Expandir metadatos de SRA

Para estudios de SRA (Sequence Read Archive), es útil expandir los atributos de las muestras usando `expand_sra_attributes()`. Esto convierte la información cruda de SRA en columnas más legibles.

```r
# Expandir atributos de SRA
rse_gene_SRP009615 <- expand_sra_attributes(rse_gene_SRP009615)

# Ver los atributos expandidos
colData(rse_gene_SRP009615)[
    ,
    grepl("^sra_attribute", colnames(colData(rse_gene_SRP009615)))
]

# DataFrame with 12 rows and 4 columns
#           sra_attribute.cells sra_attribute.shRNA_expression sra_attribute.source_name sra_attribute.treatment
# SRR387777                K562                             no                    SL2933               Puromycin
# SRR387778                K562             yes, targeting SRF                    SL2934  Puromycin, doxycycline
# ...                       ...                            ...                       ...                     ...
```

---

## 3. Exploración del Objeto RSE

```r
# Dimensiones
dim(rse_gene_SRP009615)

# Nombres de assays
assayNames(rse_gene_SRP009615)

# Ver primeras filas de las cuentas
head(assay(rse_gene_SRP009615, "counts"))

# Información de los genes (rowRanges)
rowRanges(rse_gene_SRP009615)

# Metadatos de las muestras (colData)
colData(rse_gene_SRP009615)
```

---

## 4. Ejercicio: Reproducir una imagen con iSEE

**Objetivo:** Utilizar `iSEE` para reproducir la siguiente imagen (un boxplot de expresión por condición experimental).

![Ejemplo de visualización con iSEE](ruta_a_la_imagen_referencia.png)

**Pistas proporcionadas:**
*   Utiliza el **dynamic feature selection** para elegir genes.
*   Utiliza información de las columnas para el **eje X**.
*   Utiliza información de las columnas para los **colores**.

### 4.1 Preparar los datos para iSEE

Antes de lanzar `iSEE`, podemos simplificar los nombres de las variables para que sean más legibles.

```r
# Cargar iSEE
library("iSEE")

# Simplificar nombres de atributos para mejor visualización
colData(rse_gene_SRP009615)$treatment <- colData(rse_gene_SRP009615)$sra_attribute.treatment
colData(rse_gene_SRP009615)$shRNA <- colData(rse_gene_SRP009615)$sra_attribute.shRNA_expression
colData(rse_gene_SRP009615)$cells <- colData(rse_gene_SRP009615)$sra_attribute.cells

# Verificar los nuevos nombres
head(colData(rse_gene_SRP009615)[, c("treatment", "shRNA", "cells")])
```

### 4.2 Lanzar iSEE y configurar la visualización

```r
# Lanzar la aplicación interactiva
if(interactive()) {
    iSEE::iSEE(rse_gene_SRP009615)
}
```

### 4.3 Pasos para reproducir la imagen

Sigue estos pasos dentro de la interfaz de `iSEE`:

1.  **Seleccionar tipo de panel:** En la esquina superior izquierda, selecciona **"Column data plot"** (gráfico de datos de columnas).

2.  **Configurar el eje X:**
    *   En el panel de configuración (engranaje ⚙️) del gráfico, busca la opción **"X-axis"**.
    *   Selecciona **"Column data"**.
    *   Elige la variable **`treatment`** (o `sra_attribute.treatment`). Esta variable contiene información sobre el tratamiento experimental (ej. "Puromycin" vs "Puromycin, doxycycline").

3.  **Configurar el eje Y:**
    *   En **"Y-axis"**, selecciona **"Feature assay"**.
    *   Aparecerá un campo para seleccionar el gen. Usa el **dynamic feature selection** (selección dinámica de características) para buscar y seleccionar un gen de interés.
    *   *Nota:* Para este ejercicio, podríamos elegir un gen que sepamos que responde al tratamiento. Si no hay un gen específico indicado, podemos seleccionar cualquiera, por ejemplo, el primer gen de la lista, o buscar uno con expresión variable.

4.  **Configurar el color:**
    *   En **"Color"**, selecciona **"Column data"**.
    *   Elige la variable **`shRNA`** (o `sra_attribute.shRNA_expression`). Esta variable indica si las células expresan un shRNA (control vs targeting SRF).

5.  **Ajustes adicionales (opcional):**
    *   Puedes cambiar el tipo de gráfico a **"Violin"** o **"Boxplot"** en las opciones de "Visual parameters".
    *   Ajusta el tamaño de los puntos si es necesario.

6.  **Reproducir la imagen:**
    *   Una vez configurado, el gráfico debería mostrar la expresión del gen seleccionado (eje Y) para cada condición de tratamiento (eje X), coloreado por el tipo de shRNA.
    *   Para guardar la imagen, haz click en el icono de **"Download"** (⬇️) en la esquina superior derecha del panel y selecciona **"PDF"**.

### 4.4 Código alternativo para generar el gráfico directamente

Si prefieres no usar la interfaz interactiva, puedes generar un gráfico similar con código R:

```r
# Crear un data frame con los datos necesarios
library("ggplot2")

# Extraer datos de expresión de un gen de ejemplo (primer gen)
gen_ejemplo <- rownames(rse_gene_SRP009615)[1]
expresion_gen <- assay(rse_gene_SRP009615, "counts")[gen_ejemplo, ]

# Crear data frame para ggplot
df_plot <- data.frame(
    expression = expresion_gen,
    treatment = colData(rse_gene_SRP009615)$treatment,
    shRNA = colData(rse_gene_SRP009615)$shRNA
)

# Crear boxplot
ggplot(df_plot, aes(x = treatment, y = expression, fill = shRNA)) +
    geom_boxplot() +
    theme_bw() +
    labs(
        title = paste("Expresión de", gen_ejemplo),
        x = "Tratamiento",
        y = "Cuentas",
        fill = "shRNA"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Guardar el gráfico
ggsave("figuras/boxplot_recount3_ejemplo.pdf", width = 8, height = 6)
```

---

## 5. Resumen de Conceptos Clave

| Concepto | Descripción |
|:---|:---|
| **`recount3`** | Recurso con >700,000 muestras de RNA-seq de humano y ratón procesadas uniformemente. |
| **`available_projects()`** | Función para explorar todos los proyectos disponibles en recount3. |
| **`create_rse()`** | Descarga los datos de un proyecto y crea un objeto `RangedSummarizedExperiment`. |
| **`compute_read_counts()`** | Transforma cuentas por nucleótido a cuentas por lectura (estándar para análisis). |
| **`expand_sra_attributes()`** | Expande los metadatos de SRA en columnas legibles. |
| **`iSEE`** | Herramienta interactiva para explorar objetos `SummarizedExperiment`. |
| **Dynamic feature selection** | En iSEE, permite buscar y seleccionar genes dinámicamente. |

```r
# Guardar información de la sesión para reproducibilidad
options(width = 120)
sessioninfo::session_info()
```


        *                      +                 .                +     .       |\      _,,,---,,_
     *       .              +         .                    *             ZZZzz /,`.-'`'    -.  ;-;;,_  
      .         .  *                       .    *            *                |,4-  ) )-,_. ,\ (  `'-'            
        .     +             *                 *                              <'---''(_/--'  `-'\_)             

---