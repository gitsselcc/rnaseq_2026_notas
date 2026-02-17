
---

# Notas de Clase: Introducción a R, GitHub y el Entorno de Desarrollo

**Curso:** RNA-Seq LCG-UNAM 2026
**Repositorio de Notas:**  https://github.com/gitsselcc/rnaseq_2026_notas

## Objetivos de la Clase de Hoy

1.  Familiarizarnos con las herramientas principales del curso: R, RStudio/Positron y GitHub.
2.  Configurar nuestro entorno de trabajo local y conectarlo con GitHub.
3.  Aprender a usar proyectos y paquetes clave como `here` y `usethis` para organizar nuestro trabajo.
4.  Explorar el uso de Inteligencia Artificial (GitHub Copilot, Gemini) como asistente de programación.

---

## 1. Introducción a las Herramientas

### 1.1 R y RStudio / Positron

*   **R:** Es un lenguaje de programación gratuito y de código abierto, muy potente para análisis estadístico y bioinformática, especialmente gracias a proyectos como **Bioconductor**.
*   **RStudio / Positron:** Son Entornos de Desarrollo Integrado (IDE) que facilitan el trabajo con R.
    *   **RStudio:** El IDE clásico y muy maduro para R.
    *   **Positron:** La nueva generación de IDE, también gratuito, más moderno y con mejor soporte para proyectos que mezclan R y Python. Está basado en VS Code. El hijo concebido con amor de VS y R jajajaj

### 1.2 GitHub y Git

*   **Git:** Sistema de control de versiones. Lleva un registro de todos los cambios en nuestro código, permitiéndonos volver a versiones anteriores y colaborar sin miedo.
*   **GitHub:** Plataforma en la nube para alojar repositorios de Git. Esencial para compartir código, colaborar en proyectos y crear portafolios públicos.
    *   **GitHub Pages:** Servicio que permite convertir un repositorio de GitHub en una página web estática (como la página del curso: `lcolladotor.github.io/rnaseq_LCG-UNAM_2026`).

### 1.3 Inteligencia Artificial (IA) como Asistente

*   **GitHub Copilot:** Asistente de IA que sugiere código en tiempo real. Disponible gratis para estudiantes a través de **GitHub Education**.
*   **Configuración en RStudio/Positron:** Se puede configurar Copilot siguiendo las guías oficiales de Posit.
*   **`chattr`:** Paquete de R que permite interactuar con GitHub Copilot (u otros LLMs) desde la consola de R con `chattr::chattr_app()`.
*   **`ellmer`:** Paquete moderno del `tidyverse` para interactuar con Grandes Modelos de Lenguaje (LLMs) como GitHub Copilot o Gemini. Se puede usar para tener una conversación en vivo.
*   **Prompts (Instrucciones) Efectivas:** Ejemplos vistos en clase:
    *   "Show me how to use the `here` package with comments in Spanish"
    *   "Muestrame como usar el paquete de `SummarizedExperiment`"
    *   "How do you download data with `recount3`"
*   **`lang`:** Paquete para traducir la ayuda de R a otros idiomas usando modelos de IA locales (vía `ollama`). Permite cambiar el idioma de la documentación (ej. `?lm` aparecerá en español).

---

## 2. Configuración del Entorno de Trabajo

### 2.1 Creación del Proyecto y Estructura de Archivos

Seguimos la filosofía de organización basada en proyectos. Creamos un proyecto para nuestras notas y dentro de él, una estructura ordenada.

```r
#  PASO 1: Crear el proyecto (Ejecutar en la consola) 

# Crea un nuevo proyecto en tu escritorio. Esto generará un archivo .Rproj
usethis::create_project("~/Desktop/rnaseq_2026_notas")

# Después de ejecutar esto, RStudio se reiniciará en el nuevo proyecto.
# Es crucial trabajar SIEMPRE desde un proyecto .Rproj para que las rutas relativas funcionen.

# PASO 2: Crear archivos de código R dentro del proyecto 

# Crea un archivo R vacío llamado '01-notas.R' dentro de la carpeta 'R/'
# Este archivo contendrá el código de la primera clase.
usethis::use_r("01-notas.R")

# Creamos otro archivo de ejemplo para visualizar datos.
usethis::use_r("02-visualizar-mtcars.R")
```

### 2.2 Código de la Clase: Visualización y Reproducibilidad

Dentro del archivo `R/02-visualizar-mtcars.R` (o `01-notas.R`), escribimos y ejecutamos el siguiente código. Las librerías se cargan al inicio y se usan a lo largo del script.

```r
# ==================================================
# Script: 02-visualizar-mtcars.R
# Propósito: Demostrar el flujo de trabajo básico:
#            cargar librerías, definir rutas, crear directorios,
#            generar una figura y guardar información de la sesión.
# ==================================================

# 1. Cargar librerías necesarias 
# 'library()' carga los paquetes para poder usar sus funciones.

# 'sessioninfo' proporciona información detallada sobre la versión de R y los paquetes cargados.
# Es fundamental para la reproducibilidad.
library("sessioninfo")

# 'here' resuelve un problema común: construir rutas a archivos de forma fiable.
# Busca la raíz del proyecto (donde está el archivo .Rproj) y construye rutas relativas a ella.
# Esto hace que tu código funcione en cualquier ordenador sin importar dónde esté guardado el proyecto.
library("here")

# 'ggplot2' es el paquete estrella para crear visualizaciones de datos de forma elegante y potente.
library("ggplot2")

# 2. Mensaje inicial 
# Un simple 'print' para confirmar que el script se está ejecutando.
print("Hola, soy [Tu Nombre]") # Reemplaza con tu nombre

# 3. Definir rutas de trabajo 
# Usamos 'here()' para construir rutas a las carpetas donde guardaremos resultados.
# 'dir_plots' será la ruta completa a la carpeta 'figuras'
# 'dir_rdata' será la ruta completa a la carpeta 'processed-data'
dir_plots <- here::here("figuras")
dir_rdata <- here::here("processed-data")

# Mostramos las rutas para verificar que son correctas.
print(paste("Los gráficos se guardarán en:", dir_plots))
print(paste("Los datos procesados se guardarán en:", dir_rdata))

# 4. Crear directorios (si es q no existen)
# 'dir.create' crea una carpeta.
# 'showWarnings = FALSE' evita que salga una advertencia si la carpeta ya existe.
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

# 5. Análisis y visualización de datos 
# Vamos a usar el conjunto de datos 'mtcars' que viene incluido con R.
# Queremos visualizar cómo varía el consumo de combustible (mpg, millas por galón)
# según el número de marchas (gear).

# Creamos un gráfico de caja (boxplot) con ggplot2.
# 1. ggplot(data = mtcars, aes(...)): Inicializa el gráfico con los datos y las variables.
#    - 'group = gear': Agrupa por número de marchas.
#    - 'y = mpg': La variable en el eje Y es el consumo.
# 2. geom_boxplot(): Añade la capa del gráfico de caja.
boxplot_gear_mpg <- ggplot(mtcars, aes(group = factor(gear), y = mpg)) +
    geom_boxplot() +
    labs(
        title = "Consumo de combustible por número de marchas",
        x = "Número de marchas (gear)",
        y = "Millas por galón (mpg)"
    ) +
    theme_bw() # Un tema más limpio para el gráfico.

# Mostramos el gráfico en la ventana de plots de RStudio.
print(boxplot_gear_mpg)

# 6. Guardar resultados 
# Guardamos el gráfico como un archivo PDF.

# 'file.path' combina el directorio y el nombre del archivo de forma segura.
# 'useDingbats = FALSE' evita problemas con algunos caracteres en el PDF.
pdf(file.path(dir_plots, "mtcars_gear_vs_mpg.pdf"),
    width = 8, height = 6, # Ajustamos el tamaño
    useDingbats = FALSE
)
# El código que produce el gráfico debe ir entre 'pdf()' y 'dev.off()'
print(boxplot_gear_mpg)
# 'dev.off()' cierra el dispositivo PDF y guarda el archivo.
dev.off()

# Mensaje de confirmación
print(paste("Gráfico guardado en:", file.path(dir_plots, "mtcars_gear_vs_mpg.pdf")))

# 7. Información de la sesión para reproducibilidad 
# Guardamos la información de todos los paquetes y versiones usadas.
# Esto es como una "receta" de los ingredientes exactos que usamos.
options(width = 120) # Para que la salida sea más ancha y legible.
sessioninfo::session_info()

# Opcional: Guardar esta información en un archivo de texto.
sink(file.path(dir_rdata, "session_info.txt"))
sessioninfo::session_info()
sink()
# O en formato RData para poder cargarlo después si es necesario.
save(sessioninfo::session_info(), file = file.path(dir_rdata, "session_info.RData"))

print("--- Fin del script ---")
```

---

## 3. Conexión con GitHub: Flujo de Trabajo

Para conectar nuestro proyecto local con un repositorio remoto en GitHub, seguimos estos pasos con el paquete `usethis`.

```r
# --- PASO 1: Configurar Git localmente (solo una vez por ordenador) ---
# Abre el archivo de configuración global de Git (.gitconfig)
usethis::edit_git_config()
# En el archivo que se abre, añade tu nombre y email (el asociado a tu cuenta de GitHub)
# [user]
#     name = Tu Nombre Completo
#     email = tu.email@ejemplo.com

# --- PASO 2: Generar y guardar un Token de Autenticación (PAT) ---
# Esto es necesario por razones de seguridad. El token es como una contraseña especial.

# 2a. Generar el token en el navegador (se abrirá una página de GitHub)
usethis::create_github_token()
#   - Dale un nombre descriptivo al token (ej. "LCG-UNAM-Notas").
#   - Selecciona los permisos (los que vienen por defecto suelen ser suficientes).
#   - Copia el token generado (una cadena de 40 caracteres) AHORA MISMO, solo se ve una vez.

# 2b. Guardar el token en tu ordenador de forma segura.
# Para Mac/Windows:
gitcreds::gitcreds_set()
#   - Pega el token cuando te lo pida.

# Para Linux (como en el servidor de la LCG-UNAM):
#   - Abre el archivo .Renviron
usethis::edit_r_environ()
#   - Añade una línea con tu token. Ejemplo:
#     GITHUB_PAT=ghp_tu_token_de_40_caracteres_copiado
#   - Guarda el archivo y reinicia la sesión de R (Session -> Restart R).

# --- PASO 3: Inicializar Git y conectarlo con GitHub ---

# 3a. Inicializar el repositorio Git en el proyecto local.
# Esto crea la carpeta oculta '.git' y empieza a trackear los cambios.
usethis::use_git()
#   - Te preguntará si quieres hacer un 'commit' inicial. ¡Selecciona "Sí"!

# 3b. Conectar el repositorio local con un repositorio remoto en GitHub.
# Esto crea un nuevo repositorio en tu cuenta de GitHub con el mismo nombre que tu proyecto.
usethis::use_github()
#   - Te pedirá confirmar el nombre del repositorio. Presiona "Enter" para aceptar el que te sugiere.

# ¡Listo! Tu proyecto está ahora en GitHub. La URL se mostrará en la consola.
# Ejemplo: https://github.com/tu_usuario/rnaseq_2026_notas

# A partir de ahora, puedes hacer 'commit' de tus cambios locales y 'push' para subirlos a GitHub.
```

## Resumen de Conceptos Clave

*   **`usethis`:** Un paquete increíblemente útil para automatizar tareas repetitivas en la creación y configuración de proyectos R. Lo usamos para crear archivos (`use_r`), configurar Git (`use_git`, `use_github`), editar archivos de configuración (`edit_r_environ`, `edit_git_config`) y mucho más.
*   **`here`:** La solución al problema de las rutas de archivos. Olvídate de `setwd()`. Al trabajar en un proyecto (`.Rproj`), `here::here()` siempre encuentra la ruta a la raíz de tu proyecto, permitiéndote construir rutas relativas de forma fiable (`here::here("figuras", "mi_grafico.pdf")`).
*   **`sessioninfo`:** La herramienta para la reproducibilidad. `sessioninfo::session_info()` guarda una instantánea de tu entorno de trabajo (versión de R, sistema operativo, y todos los paquetes cargados con sus versiones). Esto permite que cualquiera (o tú mismo en el futuro) pueda replicar el análisis exactamente.
*   **Proyectos (.Rproj):** La base de un flujo de trabajo ordenado. Un proyecto encapsula todo tu trabajo: código, datos, resultados y configuraciones. Es el punto de partida para cualquier análisis serio y bien planificado.

---

        *                      +                 .                +     .       |\      _,,,---,,_
     *       .              +         .                    *             ZZZzz /,`.-'`'    -.  ;-;;,_  
      .         .  *                       .    *            *                |,4-  ) )-,_. ,\ (  `'-'            
        .     +             *                 *                              <'---''(_/--'  `-'\_)             

