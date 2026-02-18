
---

# Notas de Clase: Objetos de Bioconductor para Datos de Expresión

**Curso:** RNA-Seq LCG-UNAM 2026

## Objetivos de la Clase de Hoy

1.  Comprender la estructura y utilidad del objeto `SummarizedExperiment`, la clase fundamental de Bioconductor para datos de expresión.
2.  Aprender a crear y manipular objetos `SummarizedExperiment`.
3.  Introducir `GenomicRanges` para trabajar con coordenadas genómicas.
4.  Explorar datos de forma interactiva con el paquete `iSEE`.
5.  Aplicar los conceptos con datos reales de transcriptómica espacial usando `spatialLIBD`.

---

## 1. Fundamentos: `SummarizedExperiment`

`SummarizedExperiment` es la estructura de datos por excelencia en Bioconductor para almacenar y manipular datos de experimentos de alto rendimiento (como RNA-seq o ChIP-seq). Su diseño sigue la figura 2 del artículo de referencia de Bioconductor 2015:

*   **`assays`**: Una lista de matrices 2D (filas = genes/features, columnas = muestras). Puede contener múltiples capas: cuentas crudas, CPM, log-CPM, etc.
*   **`rowRanges`**: Información genómica de las filas (ej. cromosoma, inicio, fin, strand). Usa el objeto `GRanges` del paquete `GenomicRanges`.
*   **`colData`**: Metadatos de las columnas (condición experimental, sexo, batch, etc.). Es un `DataFrame`.
*   **`metadata`**: Metadatos adicionales del experimento a nivel global (opcional).

Esta estructura garantiza que todos los datos relacionados (las cuentas, la anotación de genes y la información de las muestras) viajen juntos, evitando errores y facilitando la reproducibilidad.

```r
# Paquetes necesarios
library("SummarizedExperiment")
library("GenomicRanges")
```

### 1.1 Creación de un Objeto `SummarizedExperiment` desde Cero

Construimos un objeto con datos simulados para 200 genes y 6 muestras.

```r
# --- 1. Definir dimensiones y fijar semilla para reproducibilidad ---
nrows <- 200  # Número de genes (filas)
ncols <- 6    # Número de muestras (columnas)
set.seed(20210223)

# --- 2. Crear la matriz de cuentas (assays) ---
# Generamos una matriz de números aleatorios que simulan cuentas de expresión
counts <- matrix(runif(nrows * ncols, 1, 1e4), nrows)
rownames(counts) <- paste0("gene_", 1:nrows)
colnames(counts) <- LETTERS[1:ncols]

# --- 3. Crear la información genómica de los genes (rowRanges) ---
# Usamos GRanges para especificar cromosoma, posición y hebra
rowRanges <- GRanges(
    # Los primeros 50 genes en chr1, los siguientes 150 en chr2
    seqnames = rep(c("chr1", "chr2"), c(50, 150)),
    # Rangos IRanges: posición de inicio aleatoria y ancho fijo de 100 pb
    ranges = IRanges(floor(runif(200, 1e5, 1e6)), width = 100),
    # Hebra: '+' o '-' aleatoria
    strand = sample(c("+", "-"), 200, TRUE),
    # Metadatos adicionales por gen (feature_id)
    feature_id = sprintf("ID%03d", 1:200)
)
# Asignamos nombres a los genes para poder identificarlos
names(rowRanges) <- paste0("gene_", seq_len(length(rowRanges)))

# --- 4. Crear los metadatos de las muestras (colData) ---
# Creamos un DataFrame con la condición experimental de cada muestra
colData <- DataFrame(
    Treatment = rep(c("ChIP", "Input"), 3),
    # Asignamos nombres de fila que coinciden con las columnas de la matriz 'counts'
    row.names = LETTERS[1:6]
)

# --- 5. Ensamblar todo en un objeto SummarizedExperiment ---
# Usamos SimpleList para el assay, que es un contenedor estándar en Bioconductor
rse <- SummarizedExperiment(
    assays = SimpleList(counts = counts),
    rowRanges = rowRanges,
    colData = colData
)

# --- 6. Explorar el objeto resultante ---
rse
# class: RangedSummarizedExperiment 
# dim: 200 6 
# assays(1): counts
# rowRanges: GRanges con metadata 'feature_id'
# colnames(6): A B C D E F
# colData names(1): Treatment

# Número de genes y muestras
dim(rse)
# [1] 200   6

# Nombres de las filas (genes) y columnas (muestras)
dimnames(rse)

# Nombres de los assays disponibles
assayNames(rse)
# [1] "counts"

# Ver las primeras filas de la matriz de cuentas
head(assay(rse, "counts"))

# Acceder a la información genómica completa de los genes
rowRanges(rse)
# GRanges object with 200 ranges and 1 metadata column

# Acceder solo a los metadatos de los genes (DataFrame plano)
rowData(rse)  # Equivalente a mcols(rowRanges(rse))

# Acceder a los metadatos de las muestras
colData(rse)
# DataFrame with 6 rows and 1 column
#     Treatment
#   <character>
# A        ChIP
# B       Input
# C        ChIP
# D       Input
# E        ChIP
# F       Input
```

---

## 2. Ejercicio: Manipulación de `SummarizedExperiment`

**Instrucción:** Explica qué sucede en los siguientes comandos de R.

```r
# Comando 1
rse[1:2, ]
```

**Explicación:** Este comando utiliza el sub-setting estándar de R ( `[filas, columnas]` ) sobre un objeto `SummarizedExperiment`. Al indicar `1:2` en la posición de filas, estamos seleccionando **los primeros dos genes**. Al dejar el espacio de las columnas vacío (` , ` ), estamos seleccionando **todas las muestras**. El resultado es un nuevo objeto `RangedSummarizedExperiment` que contiene solo los datos de los genes `gene_1` y `gene_2` para las 6 muestras.

```r
# Comando 2
rse[, c("A", "D", "F")]
```

**Explicación:** Aquí se invierte la selección. El espacio de filas está vacío (`[, ...]`), lo que significa que se seleccionan **todos los genes**. En la posición de columnas, se seleccionan por nombre las muestras **"A", "D" y "F"**. El resultado es un nuevo objeto con los 200 genes, pero únicamente para esas 3 muestras específicas.

**Conclusión:** La herencia de `SummarizedExperiment` permite usar la sintaxis familiar de matrices para filtrar objetos complejos de forma intuitiva.

---

## 3. Exploración Interactiva con `iSEE`

`iSEE` (Interactive SummarizedExperiment Explorer) es un paquete de Bioconductor que proporciona una interfaz gráfica interactiva (shiny) para explorar objetos `SummarizedExperiment` y `SingleCellExperiment` sin escribir código.

```r
# Cargar el paquete
library("iSEE")

# Lanzar la aplicación con nuestro objeto 'rse'
if(interactive()){
    iSEE::iSEE(rse)
}
```

**Funcionalidades principales:**
*   **Heatmap:** Visualiza patrones de expresión.
*   **Reduced dimension plot (PCA, t-SNE, UMAP):** Explora la relación entre muestras.
*   **Column data plot:** Crea gráficos de barras, boxplots, etc., basados en los metadatos de las muestras.
*   **Feature assay plot:** Grafica la expresión de uno o varios genes.
*   **Row statistics table:** Tabla con estadísticas de los genes.

---

## 4. Ejercicio Práctico con Datos Reales: `spatialLIBD`

Vamos a trabajar con datos de transcriptómica espacial del cerebro humano (DLPFC) disponibles en el paquete `spatialLIBD`. Estos datos están almacenados en un objeto `SingleCellExperiment`, una extensión de `SummarizedExperiment` para datos de célula única.

```r
# Instalar paquete si es necesario
# BiocManager::install("spatialLIBD")

# Cargar librerías necesarias
library("spatialLIBD")
library("lobstr")  # Para ver el tamaño del objeto
library("iSEE")

# --- Descargar los datos ---
# fetch_data() descarga y cachea los datos automáticamente
sce_layer <- spatialLIBD::fetch_data("sce_layer")
# adding rname 'https://www.dropbox.com/...'

# --- Explorar el objeto descargado ---
sce_layer
# class: SingleCellExperiment 
# dim: 22331 76 
# assays(2): counts logcounts
# rownames(22331): ENSG00000243485 ENSG00000238009 ...
# rowData names(10): source type ...
# colnames(76): 151507_Layer1 151507_Layer2 ... 151676_Layer6 151676_WM
# colData names(13): sample_name layer_guess ...
# reducedDimNames(6): PCA TSNE_perplexity5 ...

# --- Verificar el tamaño en memoria ---
lobstr::obj_size(sce_layer)
# 33.99 MB  (Es un objeto de tamaño moderado)

# --- Lanzar iSEE para exploración interactiva ---
if(interactive()){
    iSEE::iSEE(sce_layer)
}
```

### 4.1 Resolución de las Preguntas del Ejercicio

Utilizando la aplicación interactiva `iSEE` con los datos `sce_layer`, exploramos los genes solicitados.

**Pregunta 1:** Explora con un heatmap la expresión de los genes *MOBP*, *MBP* y *PCP4*. Si hacemos un clustering, ¿cuáles genes se parecen más?

*   **Estrategia en iSEE:**
    1.  En el panel de `RowData`, buscar los genes por su símbolo o Ensembl ID:
        *   *MOBP* = `ENSG00000168314`
        *   *PCP4* = `ENSG00000183036`
        *   *MBP* = `ENSG00000197971`
    2.  Seleccionar un panel de tipo **"Heatmap"**.
    3.  En "Row data", elegir la opción "Selected" para usar solo los genes seleccionados.
    4.  Habilitar el clustering tanto para filas como para columnas.

*   **Resultado observado:**
    Al visualizar el heatmap y activar el clustering jerárquico, se observa que **MOBP y MBP** tienen un patrón de expresión muy similar y se agrupan en la misma rama del dendrograma. *PCP4* muestra un patrón claramente distinto y se separa de los otros dos.

*   **Interpretación biológica:** *MOBP* (Myelin-Associated Oligodendrocyte Basic Protein) y *MBP* (Myelin Basic Protein) son genes clásicos de oligodendrocitos y están involucrados en la formación de mielina. Es lógico que tengan patrones de expresión espacial similares en el cerebro. *PCP4* (Purkinje Cell Protein 4) es un gen asociado a neuronas, con una función y localización diferente.

**Pregunta 2:** ¿En qué capas se expresan más los genes *MOBP* y *MBP*?

*   **Estrategia en iSEE:**
    1.  Seleccionar un panel de tipo **"Column data plot"** (o un Feature assay plot).
    2.  Configurar el eje X para que muestre la variable de interés, por ejemplo, `layer_guess` o `layer_guess_reordered_short`, que contienen la identidad de la capa cortical (L1 a L6) y la materia blanca (WM).
    3.  Configurar el eje Y para mostrar la expresión del gen `ENSG00000168314` (MOBP) o `ENSG00000197971` (MBP). Usar `assay = "logcounts"` para ver los datos normalizados.
    4.  El tipo de gráfico puede ser un boxplot o violin plot.

*   **Resultado observado:**
    Los boxplots revelan que tanto **MOBP como MBP alcanzan su máxima expresión en la materia blanca (WM)** y en las capas más profundas de la corteza, particularmente la **capa 6 (Layer6)**. Esta es una observación biológicamente coherente, ya que la materia blanca está compuesta principalmente por axones mielinizados (producidos por oligodendrocitos) y la capa 6 es la capa cortical más interna, rica en fibras de proyección.

```r
# Código alternativo para verificar en consola (menos interactivo pero reproducible)
library("scater") # Para funciones de plotting de SCE

# Añadir los logcounts si no existen
sce_layer <- logNormCounts(sce_layer)

# Crear boxplots para MOBP
plotExpression(sce_layer, 
               features = "ENSG00000168314", 
               x = "layer_guess_reordered_short",
               colour_by = "layer_guess_reordered_short") +
    ggtitle("Expresión de MOBP por capa cortical")

# Crear boxplots para MBP
plotExpression(sce_layer, 
               features = "ENSG00000197971", 
               x = "layer_guess_reordered_short",
               colour_by = "layer_guess_reordered_short") +
    ggtitle("Expresión de MBP por capa cortical")
```

---

## 5. Resumen de Conceptos Clave

| Concepto | Descripción |
|:---|:---|
| **`SummarizedExperiment`** | Objeto contenedor que integra `assays` (datos), `rowRanges` (anotación genómica) y `colData` (metadatos de muestras). |
| **`RangedSummarizedExperiment`** | Subtipo de `SummarizedExperiment` donde `rowRanges` es obligatorio y almacena coordenadas genómicas. Es el más común. |
| **`GRanges`** | Objeto del paquete `GenomicRanges` para almacenar y manipular intervalos genómicos (cromosoma, inicio, fin, strand). |
| **`SingleCellExperiment`** | Extensión de `SummarizedExperiment` para datos de célula única. Incluye `reducedDims` para almacenar reducciones de dimensionalidad (PCA, t-SNE, UMAP). |
| **`iSEE`** | Aplicación interactiva para explorar objetos `SummarizedExperiment` y `SingleCellExperiment` sin programar. |
| **`spatialLIBD`** | Paquete con datos de transcriptómica espacial del cerebro humano. Útil para aprender y explorar. |

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