
---

# Notas de Clase: Introducción a Bioconductor

**Curso:** RNA-Seq LCG-UNAM 2026
**Fecha de la Clase:** 

## Objetivos de la Clase de Hoy

1.  Entender qué es Bioconductor y por qué es fundamental en bioinformática.
2.  Aprender a navegar por la estructura de Bioconductor: paquetes, ramas `release` y `devel`.
3.  Identificar los diferentes tipos de paquetes: Software, Annotation, Experiment Data, Workflows.
4.  Explorar cómo encontrar paquetes relevantes usando `biocViews`.
5.  Realizar un ejercicio práctico de exploración de nuevos paquetes en Bioconductor 3.22.

---

## 1. ¿Qué es Bioconductor?

Bioconductor es un repositorio de paquetes de R de código abierto y libre, especializado en el **análisis y comprensión de datos genómicos de alto rendimiento**. Es, junto con CRAN, uno de los repositorios más importantes para R.

*   **Sitio web oficial:** [http://bioconductor.org/](http://bioconductor.org/)
*   **Filosofía:** Proveer herramientas estadísticas y visuales robustas, con documentación exhaustiva (vignettes) y un fuerte enfoque en reproducibilidad.
*   **Gobernanza:** Cuenta con un *Community Advisory Board* del cual el instructor (Leo) fue parte (2020-2023).

### Artículos y Libros Fundamentales

| Año | Tipo | Título/Descripción |
|:---:|:---|:---|
| 2004 | Artículo | "Bioconductor: Open software development for computational biology and bioinformatics" (Genome Biology) |
| 2005 | Libro | "Bioinformatics and Computational Biology Solutions Using R and Bioconductor" (Springer) |
| 2015 | Artículo | "Orchestrating high-throughput genomic analysis with Bioconductor" (Nature Methods) |
---

## 2. Estructura y Navegación en Bioconductor

### 2.1 Tipos de Paquetes

Bioconductor organiza sus paquetes en cuatro categorías principales:

| Tipo de Paquete | Descripción | Ejemplo de URL |
|:---|:---|:---|
| **Software** | Paquetes que implementan métodos de análisis, algoritmos, visualizaciones, etc. Es el tipo más común. | [bioc Software](http://bioconductor.org/packages/release/bioc/) |
| **Annotation** | Paquetes que contienen datos de anotación genómica (ej. identificadores de genes, posiciones cromosómicas, rutas metabólicas) para facilitar la interpretación biológica. | [bioc Annotation](http://bioconductor.org/packages/release/data/annotation/) |
| **Experiment Data** | Paquetes que contienen conjuntos de datos de experimentos reales. Útiles para reproducir análisis de artículos o probar flujos de trabajo. | [bioc Experiment Data](http://bioconductor.org/packages/release/data/experiment/) |
| **Workflows** | Paquetes que documentan y demuestran, a través de *vignettes*, cómo combinar múltiples paquetes de Bioconductor para realizar un análisis completo (ej. RNA-seq, ChIP-seq). | [bioc Workflows](http://bioconductor.org/packages/release/workflows/) |

### 2.2 Descubriendo Paquetes con `biocViews`

`biocViews` es un sistema de clasificación jerárquica (como un árbol de categorías) que permite encontrar paquetes por temática.

*   **Árboles principales:** `Software`, `Annotation`, `ExperimentData`, `Workflow`.
*   **Ramas temáticas:** Dentro de cada árbol, los paquetes se etiquetan con términos como `Visualization`, `RNASeq`, `DifferentialExpression`, etc. Un paquete puede tener múltiples etiquetas.
*   **Ejemplo de navegación:** Software → `WorkflowStep` → `Visualization` → [Ver paquetes de visualización](http://bioconductor.org/packages/release/BiocViews.html#___Visualization)

### 2.3 Crecimiento de Bioconductor

El número de paquetes de software crece constantemente con cada lanzamiento semestral.

| Versión BioC | Período | # Paquetes Software |
|:---:|:---|:---:|
| 3.11 | Abr-Oct 2020 | 486 |
| 3.12 | Oct 2020-Abr 2021 | 506 |
| 3.14 | Oct 2021-Abr 2022 | 536 |
| 3.16 | Oct 2022-Abr 2023 | 542 |
| 3.18 | Oct 2023-Abr 2024 | 558 |
| 3.20 | Oct 2024-Abr 2025 | 566 |
| **3.22** | **Oct 2025-Abr 2026** | **588** |

### 2.4 Anatomía de la Página de un Paquete

Cada paquete tiene su propia página web con información crucial. Ejemplo: [`recount3`](https://bioconductor.org/packages/recount3).

1.  **URL:** `https://bioconductor.org/packages/<nombre_paquete>`
2.  **Badges (Etiquetas):** Estado de compilación (OK, ERROR, WARN) en diferentes sistemas operativos (Linux, macOS, Windows). ¡Revisarlas es clave para saber si el paquete funciona!
3.  **Descripción:** Párrafo resumiendo la funcionalidad del paquete.
4.  **Citación:** Cómo citar el paquete si lo usas en una publicación.
5.  **Instalación:** Comando para instalar el paquete (ver sección 3).
6.  **Documentación:**
    *   **Vignettes:** Son la documentación principal. Son documentos (PDF/HTML) que explican, con ejemplos, cómo usar el paquete. ¡Siempre empieza por aquí!
    *   **Referencia de funciones:** Lista de todas las funciones del paquete.
7.  **Detalles:**
    *   **`biocViews`:** Las etiquetas de clasificación del paquete.
    *   **Dependencias:** Qué otros paquetes necesita (`Depends`, `Imports`).
    *   **URL:** Enlace al repositorio de desarrollo (ej. GitHub).
    *   **BugReports:** Dónde reportar problemas.

### 2.5 Las Ramas: `release` y `devel`

Bioconductor se actualiza cada 6 meses (abril y octubre), creando dos ramas activas simultáneamente:

| Rama | Versión Actual (a feb 2026) | Descripción |
|:---|:---:|:---|
| **`release`** | **3.22** | Versión estable. Recomendada para la mayoría de usuarios. Es la que se debe usar para análisis publicables. |
| **`devel`** | **3.23** | Versión de desarrollo. Contiene paquetes nuevos y cambios que se probarán para el próximo lanzamiento. Para desarrolladores y usuarios que quieran probar las últimas novedades. |

Puedes ver la página de un paquete en la rama `devel` cambiando la URL: `http://bioconductor.org/packages/devel/bioc/html/recount3.html`

---

## 3. Instalación de Paquetes de Bioconductor

La forma estándar de instalar paquetes de Bioconductor es usando el paquete `BiocManager`.

```r
# --- PASO 1: Instalar BiocManager (solo una vez) ---
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# --- PASO 2: Instalar paquetes de Bioconductor ---
# La función 'install()' de BiocManager instala desde CRAN, Bioc y GitHub.
BiocManager::install("<nombre_del_paquete>")

# Ejemplo: Instalar paquetes que vimos en clase
BiocManager::install("ComplexHeatmap")   # Visualización avanzada de heatmaps
BiocManager::install("SummarizedExperiment") # Estructura de datos fundamental
BiocManager::install("recount3")         # Acceso a datos públicos de RNA-seq

# Para instalar múltiples paquetes a la vez
BiocManager::install(c("ComplexHeatmap", "SummarizedExperiment", "iModMix", "cellmig"))

# --- PASO 3: Verificar la versión de Bioconductor ---
BiocManager::version()
# [1] '3.22'  # Debería mostrar la versión release

# Para actualizar todos los paquetes a las últimas versiones
BiocManager::install()
```

---

## 4. Comunidad y Soporte

Bioconductor tiene una comunidad muy activa y varios canales de ayuda:

*   **Zulip Chat:** [https://community-bioc.zulipchat.com/](https://community-bioc.zulipchat.com/) - Para discusiones en tiempo real.
*   **Sitio de Soporte (Support Site):** [https://support.bioconductor.org/](https://support.bioconductor.org/)
    *   **¡Usa las etiquetas (`tags`) adecuadas!** Ejemplo: si tienes una pregunta sobre `recount3`, usa la etiqueta `recount3`. Esto asegura que los autores del paquete reciban una notificación automática.
    *   Es un excelente lugar para *aprender* viendo cómo otros resuelven problemas.
*   **Eventos y Cursos:**
    *   [Calendario de eventos](http://bioconductor.org/help/events/)
    *   [Materiales de cursos](http://bioconductor.org/help/course-materials/)
    *   **Conferencia anual:** [BioC2026](https://bioc2026.bioconductor.org/)
    *   **Eventos en México:** [Comunidad Bioinfo (CDSB)](https://comunidadbioinfo.github.io/es/#events)
*   **Bluesky:** [@bioconductor.bsky.social](https://bsky.app/profile/bioconductor.bsky.social)

---

## 5. Ejercicio Grupal: Exploración de Nuevos Paquetes en BioC 3.22

**Instrucción:** Exploramos la lista de nuevos paquetes en el lanzamiento de Bioconductor 3.22 ([anuncio oficial](http://bioconductor.org/news/bioc_3_22_release/)). Cada persona eligió un paquete y reportó sus hallazgos.

### Paquete 1: `iModMix` (Explorado por Yael)

*   **¿Qué hace?**
    `iModMix` implementa una estrategia basada en **redes para integrar datos multi-ómicos** (transcriptómica, proteómica, metabolómica) obtenidos de las mismas muestras.

*   **¿Cómo funciona?**
    1.  **Análisis por capa:** Primero analiza cada conjunto de datos por separado, infiriendo redes de correlación parcial usando modelos gráficos gaussianos con *graphical lasso*. Esto ayuda a priorizar asociaciones directas entre variables.
    2.  **Definición de módulos:** Estas redes se agrupan usando superposición topológica y clustering jerárquico para definir **módulos** de características (ej. genes, proteínas) que están fuertemente conectadas.
    3.  **Resumen de módulos:** Cada módulo se resume en una "característica propia" (*eigenfeature*), que es el primer componente principal de las variables del módulo.
    4.  **Integración entre capas:** Finalmente, correlaciona estas características propias entre los diferentes conjuntos de datos (ej. correlación entre un módulo de transcriptómica y uno de proteómica) para generar una **red de módulo a módulo**.

*   **Utilidad:** Permite descubrir programas biológicos coordinados a través de múltiples capas moleculares, incluso sin anotaciones funcionales previas. Es una herramienta exploratoria pero con base estadística sólida.

### Paquete 2: `cellmig` (Explorado por Marina)

*   **¿Qué hace?**
    `cellmig` es un paquete que implementa **modelos jerárquicos bayesianos** para analizar datos de **migración celular** provenientes de experimentos de imágenes de alto rendimiento (*high-throughput imaging*).

*   **¿Cuál es su innovación?**
    Maneja explícitamente la estructura anidada de los datos (células individuales dentro de réplicas técnicas, dentro de réplicas biológicas) y el ruido técnico asociado. Esto permite cuantificar con precisión los cambios en la velocidad de migración inducidos por condiciones experimentales (ej. fármacos), superando limitaciones de métodos tradicionales que ignoran estas dependencias.

*   **Características adicionales:**
    *   Compatible con R 4.5 en Linux, Windows y macOS.
    *   Permite generar **datos sintéticos**, lo cual es muy útil para optimizar el diseño experimental antes de realizar el experimento real.

### Paquete 3: `ComplexHeatmap` (Explorado por Gissel, osea yo :D)

*   **¿Qué hace?**
    `ComplexHeatmap` (versión 2.26.1) es un paquete para crear **heatmaps complejos y altamente personalizables**. Su objetivo es visualizar asociaciones entre múltiples fuentes de datos y descubrir patrones.

*   **Características principales:**
    *   Organizar múltiples heatmaps de forma flexible.
    *   Incorporar una gran variedad de gráficos de anotación: barplots, boxplots, histogramas, imágenes, textos, etc.
    *   Integrar dendrogramas, títulos y anotaciones en un solo objeto visual.

*   **Estructura orientada a objetos:**
    *   `Heatmap-class`: Representa un heatmap individual.
    *   `HeatmapList-class`: Una lista o colección de heatmaps y anotaciones.
    *   `HeatmapAnnotation-class`: Define anotaciones para filas o columnas.
    *   `ColorMapping-class`: Controla cómo se mapean los valores a colores.

*   **Autor/Mantenedor:** Zuguang Gu.
*   **Licencia:** MIT.

### Discusión en Equipo: Conclusiones

*   **Diversidad de aplicaciones:** Los paquetes reflejan la amplitud de la bioinformática moderna: desde integración multi-ómica (`iModMix`) hasta análisis especializados de imágenes celulares (`cellmig`) y visualización de datos (`ComplexHeatmap`).
*   **Enfoque estadístico:** Hay una tendencia clara hacia métodos estadísticos robustos (redes gráficas, modelos bayesianos jerárquicos) que modelan explícitamente la complejidad de los datos biológicos (estructuras anidadas, ruido técnico).
*   **Documentación:** Todos los paquetes tienen una documentacion muy buena pero es evidente que `ComplexHeatmap` es quien tiene la docuemntacion mas extensa y detallada.
*   **Conclusion:** Estos paquetes y muchos mas dentro de *BioC* muy probablemente nos ayuden a visualizar y manejar datos a gran escala.

---

## Resumen de Conceptos Clave

| Concepto | Descripción |
|:---|:---|
| **Bioconductor** | Repositorio de paquetes R para análisis de datos genómicos de alto rendimiento. |
| **Tipos de paquetes** | Software (métodos), Annotation (datos de anotación), Experiment Data (datos de ejemplo), Workflows (flujos completos). |
| **`biocViews`** | Sistema de etiquetas jerárquico para encontrar paquetes por tema (ej. `RNASeq`, `Visualization`). |
| **`release` vs `devel`** | `release` (estable, actual 3.22) para análisis; `devel` (en desarrollo, 3.23) para probar novedades. |
| **Página del paquete** | Fuente central de información: badges de estado, descripción, citación, instalación, y lo más importante: **vignettes**. |
| **Vignettes** | Documentación extensa con ejemplos de uso. ¡Siempre el primer lugar para aprender a usar un paquete! |
| **`BiocManager`** | Paquete esencial para instalar y actualizar paquetes de Bioconductor. |
| **Comunidad** | Soporte vía Zulip, support.bioconductor.org (¡usar etiquetas!!!!!)conferencias y cursos. |

```r
# Código de ejemplo: Verificar nuestra versión de Bioconductor
library(BiocManager)
version()
```

        *                      +                 .                +     .       |\      _,,,---,,_
     *       .              +         .                    *             ZZZzz /,`.-'`'    -.  ;-;;,_  
      .         .  *                       .    *            *                |,4-  ) )-,_. ,\ (  `'-'            
        .     +             *                 *                              <'---''(_/--'  `-'\_)             

---