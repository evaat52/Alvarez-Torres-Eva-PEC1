# 1. Descargar librerías necesarias

library(httr) # Para descargar los archivos
library(readr)  # Para leer archivos CSV
library(tidyverse) # Para manipulación de datos
library(S4Vectors) #Para funcion DataFrame
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("BiocManager") # Gestor de paquetes de Bioconductor (instalación de SummarizedExperiment)
# Instalar el paquete SummarizedExperiment
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment) # manejar objetos de datos con metadatos (bioconductor)
library(reshape2)  # Para reestructurar data frames (formato largo ↔ ancho)
library(ggpubr). #gráficos estadísticos
library(ggplot2)  #gráficos estadísticos
library(scales) #personalizar escalas de ejes en ggplot2
library(pheatmap) #mapa de calor
library(dplyr) #manipulación  de data frames

# 2. Descargar el archivo, ver su estructura y las primeras filas

# URL del archivo CSV en GitHub en versión raw
file_url <- "https://raw.githubusercontent.com/nutrimetabolomics/metaboData/refs/heads/main/Datasets/2018-MetabotypingPaper/DataValues_S013.csv"
# Descargar el archivo CSV
download.file(file_url, destfile = "metabolomics_dataset.csv", mode = "wb")
# Leer el archivo CSV en R
metabolomics_data <- read.csv("metabolomics_dataset.csv")
#Estructura del dataset
str(metabolomics_data)
# Ver las primeras filas del dataset
head(metabolomics_data)
# Verifica el número de filas y columnas
dim(metabolomics_data) 

# 3. Creación del SummarizedExperiment Object

# Extraer la información de los pacientes
colData <- metabolomics_data[, 1:6]  # Variables cualitativas
colData <- DataFrame(colData)  # Convertir a DataFrame de Bioconductor
dim(colData)  
head(colData)

# Extraer la matriz de datos (variables cuantitativas)
assay_data <- as.matrix(metabolomics_data[, 7:ncol(metabolomics_data)])
dim(assay_data)  
head(assay_data[, 1:5])  # Para ver los primeros valores

# Transponer la matriz de datos para que las muestras estén en filas
assay_data <- t(assay_data)

# Verificar dimensiones después de la transposición
dim(assay_data)  

#  SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(counts = assay_data),
  colData = colData
)
se

# Ver la estructura del objeto SummarizedExperiment
str(se)

# Ver la estructura de los metadatos del objeto SummarizedExperiment
str(colData(se))

# Ver la estructura de colData directamente
colData(se)

# 4. Análisis de los metadatos

# Extraer metadatos
meta <- as.data.frame(colData(se))

# Tablas básicas
table(meta$SURGERY)
table(meta$Group)
table(meta$GENDER)

# Histograma de edad
ggplot(meta, aes(x = as.numeric(AGE))) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "white") +
  labs(title = "Distribución de edades", x = "Edad", y = "Frecuencia") +
  theme_minimal()

# Boxplots de edad por tipo de cirugía:
ggplot(meta, aes(x = SURGERY, y = AGE, fill = SURGERY)) +
  geom_boxplot() + theme_minimal()

# Boxplots de edad por tipo de sexo:
ggplot(meta, aes(x = GENDER, y = AGE, fill = GENDER)) +
  geom_boxplot() + theme_minimal()

# 5. Estudio de homogeneidad en T0

# 5.1 Análisis de componentes principales en T0
# Extraer la matriz de datos de los metabolitos (o un subconjunto)
vars_T0 <- grep("_T0$", rownames(se), value = TRUE)
matriz_T0 <- t(assay(se)[vars_T0, ])  # Pacientes en filas, variables en columnas

# Quitar columnas con NA
matriz_T0 <- matriz_T0[, colSums(is.na(matriz_T0)) == 0]

# PCA con centrado y escalado
pca_res <- prcomp(matriz_T0, scale. = TRUE)

# Convertir a data frame para ggplot
pca_df <- as.data.frame(pca_res$x)
pca_df$Surgery <- colData(se)$SURGERY
pca_df$Group <- as.factor(colData(se)$Group)

# Plot con ggplot2
library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Surgery, shape = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA sobre metabolitos en T0",
       x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "%)"),
       y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "%)")) +
  theme_minimal()

# 5.2 Creación de mapa de calor
# Selección y normalización
matriz <- assay(se)
desviaciones <- apply(matriz, 1, sd, na.rm = TRUE)
top30 <- names(sort(desviaciones, decreasing = TRUE))[1:30]
matriz_top30 <- matriz[top30, ]
matriz_top30 <- matriz_top30[apply(matriz_top30, 1, function(x) all(!is.na(x))), ]
matriz_escalada <- t(scale(t(matriz_top30)))

#  Heatmap 
pheatmap(matriz_escalada,
         show_rownames = TRUE,
         main = "Heatmap: Top 30 metabolitos más variables",
         fontsize_row = 6)

# 5.3 Prueba t para cada metabolito en T0 entre tipos de cirugía
metabolitos_T0 <- grep("_T0$", rownames(se), value = TRUE)
resultados <- lapply(metabolitos_T0, function(met) {
  valores <- assay(se)[met, ]
  cirugia <- colData(se)$SURGERY
  t_test <- tryCatch(t.test(valores ~ cirugia), error = function(e) NULL)
  if (!is.null(t_test)) {
    data.frame(
      metabolito = met,
      p_valor = t_test$p.value
    )
  }
})
resultados_df <- do.call(rbind, resultados)
head(resultados_df[order(resultados_df$p_valor), ], 10)
