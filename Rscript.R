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