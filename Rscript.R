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
