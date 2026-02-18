# README: Curso de Introducción a RNA-seq

## Análisis de datos de RNA-seq con R/Bioconductor

**Autora:** Yeimi Gissel Contreras Cornejo, 2026

**Curso:** Introducción a RNA-seq - LCG-UNAM 2026

**Instructor:** Leonardo Collado-Torres

---

Este repositorio contiene las notas, ejercicios y proyectos desarrollados durante el curso de Introducción a RNA-seq impartido en la Licenciatura en Ciencias Genómicas (LCG-UNAM) en 2026. El curso proporciona una formación completa en el análisis de datos de RNA-seq utilizando el ecosistema de R y Bioconductor, desde los fundamentos hasta análisis avanzados de expresión diferencial.

---

##  Objetivos del Curso

- Comprender los fundamentos de R, RStudio/Positron y GitHub para investigación reproducible
- Dominar el uso de objetos de Bioconductor (`SummarizedExperiment`, `SingleCellExperiment`, `GRanges`)
- Aprender a acceder y procesar datos públicos de RNA-seq con `recount3`
- Desarrollar modelos estadísticos para expresión diferencial con `limma`-`voom` y `edgeR`
- Visualizar resultados con `ggplot2`, `pheatmap` y `ComplexHeatmap`
- Explorar datos interactivamente con `iSEE`
- Aplicar los conocimientos a casos de estudio reales (lupus, Alzheimer, datos espaciales)

---

## 📁 Estructura del Repositorio

```
rnaseq_2026_notas/
│
├── README.md                    # Este archivo
├── rnaseq_2026_notas.Rproj      # Proyecto de RStudio
├── .gitignore                    # Archivos ignorados por Git
│
├── notas/                        # Notas de clase en Markdown
│   ├── 01_IntroduccionRGithubED.md
│   ├── 02_Bioconductor.md
│   ├── 03_ObjetosBioconductor.md
│   ├── 04_Recount.md
│   ├── 05_ModelosEstadisticos.md
│   └── Platicas_rna.md
│
├── R/                            # Scripts de código ejecutable
│   ├── 01_Clase.R
│   ├── 02_Clase.R
│   ├── 03_Clase.R
│   ├── 03_SummarizedExperiment.R
│   ├── 04_Clase.R
│   ├── 04_recount3.R
│   └── 05_modelos.R
│
├── figuras/                      # Gráficos y visualizaciones generadas
│   ├── mtcars_gear_vs_mpg.pdf
│   ├── heatmap_basico.pdf
│   ├── heatmap_intermedio_anotaciones.pdf
│   ├── pheatmap_con_nombres.pdf
│   ├── ReducedDimensionPlot1.pdf
│   └── SecondFeatureAssayPlot1.pdf
│
└── processed-data/               # Datos procesados y objetos guardados
    ├── session_info.txt
    ├── session_info.RData
    ├── SummarizedExperiment_ejemplo.rds
    ├── sesion_bioconductor_complexheatmap.txt
    ├── clase2_complexheatmap.RData
    ├── rse_creado_en_clase.rds
    ├── sce_layer_spatialLIBD.rds
    ├── sesion_clase3.txt
    ├── clase3_completa.RData
    └── clase4_recount3.RData
```

---

## Contenido del Curso

### Clase 1: Introducción a R y GitHub
- Instalación de R, RStudio y Positron
- Configuración de GitHub y tokens de autenticación
- Uso de `usethis` y `here` para organización de proyectos
- Introducción a IA como asistente de programación (GitHub Copilot, `ellmer`, `lang`)
- Primer script reproducible con `sessioninfo`

### Clase 2: Introducción a Bioconductor
- Filosofía y estructura de Bioconductor
- Tipos de paquetes: Software, Annotation, Experiment Data, Workflows
- Ramas `release` y `devel`
- Exploración de paquetes con `biocViews`
- **Ejercicio:** Análisis detallado de `ComplexHeatmap` (v2.26.1)
  - Heatmaps básicos e intermedios
  - Anotaciones en columnas y filas

### Clase 3: Objetos de Bioconductor para Datos de Expresión
- `SummarizedExperiment`: estructura fundamental
- `GenomicRanges` para coordenadas genómicas
- `SingleCellExperiment` para datos de célula única
- Exploración interactiva con `iSEE`
- **Ejercicio práctico:** Datos de transcriptómica espacial con `spatialLIBD`
  - Expresión de MOBP, MBP y PCP4 en capas corticales
  - Identificación de patrones de expresión por capa

### Clase 4: Datos de RNA-seq a través de recount3
- Filosofía de democratización de datos
- Búsqueda de proyectos con `available_projects()`
- Descarga de datos con `create_rse()`
- Transformación de cuentas con `compute_read_counts()`
- Expansión de metadatos de SRA
- **Ejercicio:** Reproducción de visualizaciones con `iSEE`

### Clase 5: Modelos Estadísticos y Expresión Diferencial
- Sintaxis de fórmulas en R y `model.matrix()`
- Visualización de modelos con `ExploreModelMatrix`
- **Ejercicios resueltos:**
  - Interpretación de términos en matrices de diseño
  - Importancia del intercepto (0 en fórmulas)
- Análisis completo de expresión diferencial con `limma-voom`
  - Control de calidad y filtrado
  - Normalización con `edgeR`
  - Definición de modelo estadístico
  - Obtención de genes diferencialmente expresados
  - **Ejercicio:** Heatmap de top 50 genes con nombres (uso de `match()`)
- Visualización con MDS y volcanoplots

### Pláticas Invitadas

#### Plática 1: María Gutiérrez-Arcelus - Systems Immunology
- Genética de enfermedades inmunes (lupus)
- Variantes no codificantes y regulación génica
- Integración de GWAS, ATAC-seq y expresión
- Validación funcional del gen TNFSF4 en lupus

#### Plática 2: Daianna González Padilla - RNA-seq en el mundo real
- Fuentes de error: humanos, tecnología, metodología, biología
- Diseño experimental: teoría vs realidad
- Problemas en muestreo, extracción, preparación de librerías
- Estrategias prácticas para mitigar sesgos

#### Plática 3: Gabriel Ramírez-Vilchis - snRNA-seq del núcleo basal de Meynert
- Gen APOE y riesgo de Alzheimer (APOE4 vs APOE2)
- Vulnerabilidad del núcleo basal en Alzheimer (50% pérdida neuronal)
- snRNA-seq para estudiar composición celular y expresión
- Objeto `SingleCellSummarizedExperiment`

---

## 💻 Requisitos Técnicos

### Software
- R (>= 4.4.0)
- RStudio o Positron
- Git
- Cuenta en GitHub

### Paquetes de R
```r
# CRAN
install.packages(c(
    "here", "sessioninfo", "ggplot2", "pheatmap", 
    "RColorBrewer", "cowplot", "dplyr"
))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "SummarizedExperiment", "GenomicRanges", "SingleCellExperiment",
    "recount3", "iSEE", "ExploreModelMatrix", "spatialLIBD",
    "edgeR", "limma", "ComplexHeatmap", "circlize",
    "scater", "lobstr"
))
```

---

##  Cómo Usar Este Repositorio

### Opción 1: Clonar el repositorio
```bash
git clone https://github.com/yeimicc/rnaseq-lcg-2026.git
cd rnaseq-lcg-2026
```

### Opción 2: Descargar como ZIP
Descargar desde GitHub y descomprimir localmente.

### Opción 3: Usar `usethis`
```r
usethis::use_course("yeimicc/rnaseq-lcg-2026")
```

### Para reproducir los análisis
1. Abrir el proyecto (.Rproj) en RStudio o Positron.
2. Ejecutar los scripts en orden numérico
3. Los archivos Rmd pueden renderizarse a HTML/PDF

---

## 📊 Resultados Principales Obtenidos

### Análisis de expresión diferencial (SRP045638)
- **65 muestras** de tejido cerebral (DLPFC) tras control de calidad
- **46,929 genes** tras filtrado
- **34,079 genes diferencialmente expresados** entre prenatal y postnatal (FDR < 5%)
- Genes top: ZSCAN2, VASH2, KIAA0922

### Visualizaciones generadas
- Heatmaps con `ComplexHeatmap` y `pheatmap`
- Boxplots de calidad por grupo experimental
- MDS plots por edad y sexo
- MA-plots y volcanoplots de resultados DE

---

##  Citas y Referencias

### Paquetes principales utilizados

```r
citation("SummarizedExperiment")
citation("recount3")
citation("limma")
citation("edgeR")
citation("ComplexHeatmap")
citation("ExploreModelMatrix")
```

### Artículos clave
- Collado-Torres L, et al. (2017). "Reproducible RNA-seq analysis using recount2". *Nature Biotechnology*. DOI: 10.1038/nbt.3838
- Wilks C, et al. (2021). "recount3: summaries and queries for large-scale RNA-seq expression and splicing". *Genome Biology*. DOI: 10.1186/s13059-021-02533-6
- Law CW, et al. (2014). "voom: precision weights unlock linear model analysis tools for RNA-seq read counts". *Genome Biology*. DOI: 10.1186/gb-2014-15-2-r29
- Gu Z, et al. (2016). "Complex heatmaps reveal patterns and correlations in multidimensional genomic data". *Bioinformatics*. DOI: 10.1093/bioinformatics/btw313

---

##  Contacto

**Instructor:**
- Leonardo Collado-Torres
- Email: lcolladotor@gmail.com
- Bluesky: [@lcolladotor.bsky.social](https://bsky.app/profile/lcolladotor.bsky.social)
- GitHub: [lcolladotor](https://github.com/lcolladotor)

**Autora del repositorio:**
- Yeimi Gissel Contreras Cornejo
- Email: yeimicc@lcg.unam.mx
- GitHub: [yeimicc](https://github.com/yeimicc)

**Licenciatura en Ciencias Genómicas (LCG-UNAM)**
- https://www.lcg.unam.mx/

---

---

## Sobre Reproducibilidad

Todo el código en este repositorio incluye:
- Semillas aleatorias (`set.seed()`) para reproducibilidad
- Información de sesión (`sessioninfo::session_info()`)
- Documentación de versiones de paquetes
- Datos de acceso público (recount3, spatialLIBD)

---

##  Licencia

Este material está disponible para fines educativos. Los paquetes de R y Bioconductor utilizados tienen sus propias licencias (mayormente GPL, MIT, Artistic).

---

## 🙏 Agradecimientos

Al instructor Leonardo Collado-Torres por compartir su experiencia y materiales del curso. A los profesores invitados: María Gutiérrez-Arcelus, Daianna González Padilla y Gabriel Ramírez-Vilchis. A la LCG-UNAM por proporcionar el espacio y recursos para este curso.

---

        *                      +                 .                +     .       |\      _,,,---,,_
     *       .              +         .                    *             ZZZzz /,`.-'`'    -.  ;-;;,_  
      .         .  *                       .    *            *                |,4-  ) )-,_. ,\ (  `'-'            
        .     +             *                 *                              <'---''(_/--'  `-'\_)       
---
