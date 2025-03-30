# En este script incluyo todos los pasos realizados para la exploración de los datos, incluyendo la creación del objeto Summarizedexperiment.

# Cargo los paquetes necesarios

library(readxl)
library(S4Vectors)
library(SummarizedExperiment)

# Cargo los datos

data <- read_excel("~/OneDrive/Documentos/Master_Bioinf_Bioest/Tercer_semestre/Análisis de datos ómicos/PEC 1/DATOS REPOSITORIO/TIO2+PTYR-human-MSS+MSIvsPD.XLSX")

# Veo la estructura

str(data)

# Compruebo la presencia de valores faltantes

sum(is.na(data))

# CREACIÓN DEL OBJETO SummarizedExperiment

# Creo una matriz que solo contenga los valores de abundancia

matrix <- as.matrix(data[, -c(1, 2, 3, 4, 17, 18)])

# Asigno la secuencia como nombre de fila

rownames(matrix) <- data$SequenceModifications

# Creo metadatos con los grupos y las réplicas

colData <- data.frame(
  Sample = colnames(matrix),
  Group = rep(c("MSS", "PD"), each = 6), # Indico que cada 6 es un grupo diferente
  Replicate = rep(1:2, times = 6) # Indico que hay dos réplicas por muestra
)

# Veo la estructura

print(colData)

# Creo rowData con los datos de las demás columnas para cada fosfopéptido

rowData <- DataFrame(
  SequenceModifications = data$SequenceModifications,  # Añado también la secuencia y modificaciones (también son nombre de fila)
  Accession = data$Accession,  # Código de acceso a la proteína
  Description = data$Description,  # Descripción de la proteína
  Score = data$Score, # Score de la identificación
  Class = data$CLASS, # Clasificación
  Phosphorilation = data$PHOSPHO  # Información de fosforilación
)

# Veo la estructura de rowData

print(rowData)

# Creo el objeto SummarizedExperiment añadiendo los datos creados anteriormente

SE <- SummarizedExperiment(
  assays = list(abundance = matrix),
  colData = colData,
  rowData = rowData
)

# Veo el resumen

SE

# Guardo mi objeto SE

save(SE, file = "Objeto_SummarizedExperiment.rda")

# Guardo los datos en formato texto

write.csv(matrix, "datos.csv")