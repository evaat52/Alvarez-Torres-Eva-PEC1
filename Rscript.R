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

# Guardar el objeto 
save(se, file = "SummarizedExperiment.Rda")

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


# 6. Tendencias de metabolitos de importancia clínica y medidas antropométricas

#Tendencias de los metabolitos y medidas fundamentales de salud metabólica en función del grupo
plot_mean_se_se <- function(var_base) {
  tiempos <- c("T0", "T2", "T4", "T5")
  filas_var <- paste0(var_base, "_", tiempos)
  
  # Verificar que existen en el objeto
  filas_existentes <- filas_var[filas_var %in% rownames(se)]
  if (length(filas_existentes) < 2) {
    message("No hay suficientes tiempos para: ", var_base)
    return(NULL)
  }
  
  # Extraer datos y reorganizar
  datos <- t(assay(se)[filas_existentes, ])
  colnames(datos) <- gsub(paste0(var_base, "_"), "", filas_existentes)
  df <- cbind(as.data.frame(colData(se)), datos)
  
  # Pasar a formato largo
  df_long <- melt(df, id.vars = c("SUBJECTS", "Group", "SURGERY", "GENDER"),
                  variable.name = "Tiempo", value.name = var_base)
  # Limpiar valores no válidos en la columna Tiempo
  df_long <- df_long[df_long$Tiempo %in% tiempos, ]
  df_long$Tiempo <- factor(df_long$Tiempo, levels = tiempos)
  df_long$Group <- as.factor(df_long$Group)
  
  # Calcular medias y error estándar
  stat_df <- df_long %>%
    group_by(Tiempo, Group) %>%
    summarise(
      media = mean(.data[[var_base]], na.rm = TRUE),
      se = sd(.data[[var_base]], na.rm = TRUE)/sqrt(n())
    )
  
  # Graficar
  ggplot(stat_df, aes(x = Tiempo, y = media, group = Group, color = Group)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = media - se, ymax = media + se), width = 0.2) +
    labs(title = paste("Evolución de", var_base, "por grupo"),
         y = var_base, x = "Tiempo") +
    theme_minimal()
}

plot_mean_se_se("PESO")
plot_mean_se_se("bmi")
plot_mean_se_se("CINT")
plot_mean_se_se("CAD")
plot_mean_se_se("HBA1C")
plot_mean_se_se("LDL")
plot_mean_se_se("HDL")


#Tendencias de los metabolitos y medidas fundamentales de salud metabólica en función del tipo de cirugía
plot_mean_se_se_by_surgery <- function(var_base) {
  tiempos <- c("T0", "T2", "T4", "T5")
  filas_var <- paste0(var_base, "_", tiempos)
  
  # Verificar que existen
  filas_existentes <- filas_var[filas_var %in% rownames(se)]
  if (length(filas_existentes) < 2) {
    message("No hay suficientes tiempos para: ", var_base)
    return(NULL)
  }
  
  # Extraer datos
  datos <- t(assay(se)[filas_existentes, ])
  colnames(datos) <- gsub(paste0(var_base, "_"), "", filas_existentes)
  df <- cbind(as.data.frame(colData(se)), datos)
  
  # Pasar a formato largo
  df_long <- melt(df, id.vars = c("SUBJECTS", "Group", "SURGERY", "GENDER"),
                  variable.name = "Tiempo", value.name = var_base)
  # Limpiar valores no válidos en la columna Tiempo
  df_long <- df_long[df_long$Tiempo %in% tiempos, ]
  df_long$Tiempo <- factor(df_long$Tiempo, levels = tiempos)
  df_long$SURGERY <- as.factor(df_long$SURGERY)
  
  # Calcular medias y error estándar por cirugía
  stat_df <- df_long %>%
    group_by(Tiempo, SURGERY) %>%
    summarise(
      media = mean(.data[[var_base]], na.rm = TRUE),
      se = sd(.data[[var_base]], na.rm = TRUE)/sqrt(n())
    )
  
  # Graficar
  ggplot(stat_df, aes(x = Tiempo, y = media, group = SURGERY, color = SURGERY)) +
    geom_line(size = 1.2) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = media - se, ymax = media + se), width = 0.2) +
    labs(title = paste("Evolución de", var_base, "por tipo de cirugía"),
         y = var_base, x = "Tiempo", color = "Cirugía") +
    theme_minimal()
}
plot_mean_se_se_by_surgery("PESO")
plot_mean_se_se_by_surgery("bmi")
plot_mean_se_se_by_surgery("CINT")
plot_mean_se_se_by_surgery("CAD")
plot_mean_se_se_by_surgery("HBA1C")
plot_mean_se_se_by_surgery("LDL")
plot_mean_se_se_by_surgery("HDL")


# 6.1 Evolución de LDL según tipo de cirugía 
# Extraer valores LDL_T0 y LDL_T5
ldl_df <- data.frame(
  LDL_T0 = assay(se)["LDL_T0", ],
  LDL_T5 = assay(se)["LDL_T5", ],
  SURGERY = colData(se)$SURGERY
)

ldl_df %>%
  group_by(SURGERY) %>%
  summarise(
    media_T0 = mean(LDL_T0, na.rm = TRUE),
    media_T5 = mean(LDL_T5, na.rm = TRUE),
    diferencia = media_T5 - media_T0
  )

# Crear dataframe con LDL_T0, LDL_T5 y cirugía
ldl_test_df <- data.frame(
  LDL_T0 = assay(se)["LDL_T0", ],
  LDL_T5 = assay(se)["LDL_T5", ],
  SURGERY = colData(se)$SURGERY
)

# Calcular el cambio individual
ldl_test_df$LDL_cambio <- ldl_test_df$LDL_T5 - ldl_test_df$LDL_T0

# Realizar el t-test
t.test(LDL_cambio ~ SURGERY, data = ldl_test_df)

# ANOVA para LDL
df_ldl <- data.frame(
  LDL = c(assay(se)[paste0("LDL_T0"), ],
          assay(se)[paste0("LDL_T5"), ]),
  Tiempo = rep(c("T0", "T5"), each = ncol(se)),
  Cirugía = rep(colData(se)$SURGERY, 2)
)
anova_result <- aov(LDL ~ Tiempo * Cirugía, data = df_ldl)
summary(anova_result)

# 6.2 Metabolitos en tiempo 0 y tiempo 5 según cirugía

# Función para generar boxplots comparativos por cirugía
boxplot_cirugia <- function(var_base, tiempo) {
  var <- paste0(var_base, "_", tiempo)
  df <- data.frame(
    valor = assay(se)[var, ],
    SURGERY = colData(se)$SURGERY
  )
  
  ggplot(df, aes(x = SURGERY, y = valor, fill = SURGERY)) +
    geom_boxplot() +
    labs(title = paste("Boxplot de", var_base, "en", tiempo),
         y = var_base, x = "Cirugía") +
    theme_minimal()
}

boxplot_cirugia("LDL", "T0")
boxplot_cirugia("LDL", "T5")
boxplot_cirugia("HDL", "T0")
boxplot_cirugia("HDL", "T5")
boxplot_cirugia("HBA1C", "T0")
boxplot_cirugia("HBA1C", "T5")

# 7. Explorar si la evolución favorable de los metabolitos se correlaciona con la pérdida de peso

peso <- assay(se)[c("PESO_T0", "PESO_T5"), ]
pérdida_peso <- peso["PESO_T0", ] - peso["PESO_T5", ]
colData(se)$perdida_peso <- pérdida_peso

datos <- as.data.frame(colData(se))
cor.test(datos$perdida_peso, assay(se)["LDL_T5", ])
cor.test(datos$perdida_peso, assay(se)["HDL_T5", ])
cor.test(datos$perdida_peso, assay(se)["HBA1C_T5", ])
