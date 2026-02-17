# ==================================================
# 01-notas-clase.R - Código de la primera clase
# Curso: RNA-Seq LCG-UNAM 2026
# ==================================================

## 1. CARGAR LIBRERÍAS =====================================================
# Paquetes necesarios para el flujo de trabajo
library("sessioninfo") # Para reproducibilidad: guarda info de la sesión
library("here") # Para manejo robusto de rutas de archivos
library("ggplot2") # Para visualización de datos

## 2. CONFIGURACIÓN INICIAL ================================================
print("Hola, soy Gissel") # Mensaje de confirmación

# Definir directorios para organizar outputs
dir_plots <- here::here("figuras") # Carpeta para gráficos
dir_rdata <- here::here("processed-data") # Carpeta para datos procesados
print(paste("Gráficos en:", dir_plots))
print(paste("Datos en:", dir_rdata))

# Crear directorios si no existen
dir.create(dir_plots, showWarnings = FALSE)
dir.create(dir_rdata, showWarnings = FALSE)

## 3. ANÁLISIS Y VISUALIZACIÓN =============================================
# Usando dataset mtcars (viene con R)
# Boxplot: consumo (mpg) vs número de marchas (gear)
boxplot_gear_mpg <- ggplot(mtcars, aes(group = factor(gear), y = mpg)) +
  geom_boxplot() +
  labs(
    title = "Consumo por número de marchas",
    x = "Número de marchas (gear)",
    y = "Millas por galón (mpg)"
  ) +
  theme_bw() # Tema más limpio

# Mostrar gráfico
print(boxplot_gear_mpg)

## 4. GUARDAR RESULTADOS ===================================================
# Guardar gráfico como PDF
pdf(
  file.path(dir_plots, "mtcars_gear_vs_mpg.pdf"),
  width = 8,
  height = 6,
  useDingbats = FALSE
)
print(boxplot_gear_mpg)
dev.off() # Importante: cerrar dispositivo PDF

print(paste(
  "Gráfico guardado en:",
  file.path(dir_plots, "mtcars_gear_vs_mpg.pdf")
))

## 5. REPRODUCIBILIDAD =====================================================#

# 1. Guardar la info en un objeto
info <- sessioninfo::session_info()

# 2. Guardar ese objeto en un archivo RData
save(info, file = file.path(dir_rdata, "session_info.RData"))

# 3. Si quieres también un archivo de texto
sink(file.path(dir_rdata, "session_info.txt"))
print(info)
sink()


print("--- Script completado exitosamente ---")

## 6. CONFIGURACIÓN DE GITHUB (ejecutar una sola vez) ======================
## Estos comandos se corren en la consola, no en el script

# # Configurar usuario de Git (abre archivo para editar)
# usethis::edit_git_config()
# # Añadir en el archivo:
# # [user]
# #     name = Tu Nombre
# #     email = tu.email@ejemplo.com

# # Generar token de GitHub (abre navegador)
# usethis::create_github_token()
# # Copiar el token generado

# # Guardar token (Mac/Windows)
# gitcreds::gitcreds_set()
#
# # O para Linux (servidor LCG-UNAM):
# usethis::edit_r_environ()
# # Añadir: GITHUB_PAT=tu_token_de_40_caracteres

# # Inicializar git y conectar con GitHub
# usethis::use_git()        # Hacer commit inicial
# usethis::use_github()     # Conectar con repositorio remoto

## 7. COMANDOS ÚTILES DE IA (probados en clase) ============================
## Para usar GitHub Copilot desde R:
# chattr::chattr_app()

## Para usar ellmer con GitHub Copilot:
# library("ellmer")
# chat <- chat_github()
# live_console(chat)
# # Ejemplos: "What is my name?" "Who created R?"

## Para traducir ayuda de R a español:
# Sys.setenv(LANGUAGE = "Spanish")
# library("lang")
# llm_use("ollama", "llama3.2", seed = 100)
# ?lm  # La ayuda aparecerá en español
