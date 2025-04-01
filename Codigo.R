# En este script incluyo todos los pasos realizados para:
# 1. La creación del objeto SummarizedExperiment
# 2. La exploración y preprocesamiento de los datos
# 3. El análisis estadístico y su visualización

# Cargo los paquetes necesarios

library(readxl)
library(S4Vectors)
library(SummarizedExperiment)
library(limma)
library(edgeR)

# Cargo los datos

data <- read_excel("~/OneDrive/Documentos/Master_Bioinf_Bioest/Tercer_semestre/Análisis de datos ómicos/PEC 1/DATOS REPOSITORIO/TIO2+PTYR-human-MSS+MSIvsPD.XLSX")

# ANÁLISIS EXPLORATORIO

# Veo la estructura

str(data)

# CREACIÓN DEL OBJETO SummarizedExperiment

# Creo una matriz que contenga los valores de abundancia

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

rowData <- data.frame(
  SequenceModifications = data$SequenceModifications,  # Añado también la secuencia y modificaciones (también son nombre de fila)
  Accession = data$Accession,  # Código de acceso a la proteína
  Description = data$Description,  # Descripción de la proteína
  Score = data$Score, # Score
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

# Añado metadatos sobre el nombre del dataset, incluyendo aparte el año y el área al que pertenece

metadata(SE) <- list(
  dataset_name = "2018-Phosphoproteomics",
  year = 2018,
  research_field = "Phosphoproteomics"
)

# Veo el objeto

print(SE)

# Guardo mi objeto SE

save(SE, file = "Datos_SummarizedExperiment.rda")

# Guardo los datos en formato texto

write.csv(matrix, "datos.csv")

# ANÁLISIS EXPLORATORIO

# Compruebo la dimensión

dim(matrix)

# Veo su resumen estadístico

summary(matrix)

# Compruebo la presencia de valores faltantes

any(is.na(matrix))

# Observo la distribución de los datos en cada muestra

boxplot(t(as.matrix(assay(SE))) ~ colData$Sample,
        xlab = "Muestra", ylab = "Abundancia",
        main = "Boxplot de la abundancia de los fosfopéptidos")

# Compruebo la normalidad de los datos

qqnorm(as.vector(matrix))
qqline(as.vector(matrix), col = "brown")

# ANÁLISIS ESTADÍSTICO Y VISUALIZACIÓN

# Análisis de expresión diferencial

# Creo el objeto DGEList y normalizo los datos

dge <- DGEList(counts = assay(SE))
dge <- calcNormFactors(dge) 

# Creo el diseño de la matriz de factores (a partir de los grupos y las réplicas)

design <- model.matrix(~ Group + Replicate, data = colData)

# Aplico voom para ajustar las varianzas

v <- voom(dge, design = design, plot = TRUE)

# Realizo el análisis de expresión diferencial con limma

fit <- lmFit(v, design)

# Realizo el contraste para comparar los dos grupos (MSS vs. PD), lo aplico al modelo y ajusto (ajuste bayesiano)

contrast_matrix <- makeContrasts(GroupPD, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Obtengo los resultados

results <- topTable(fit2, adjust = "fdr", number = Inf)

# Selecciono los fosfopéptidos con p ajustados < 0.05

significant_peptides <- results[results$adj.P.Val < 0.05, ]

# Los muestro

print(significant_peptides$ID)

# Visualización de los resultados

# Creo un Volcano plot 

plot(results$logFC, -log10(results$adj.P.Val), 
     pch = 20, cex = 0.5, 
     col = ifelse(results$adj.P.Val < 0.05, "purple", "orange"), 
     xlab = "Cambio en la abundancia (Log2)", 
     ylab = "p-valor ajustado (-Log10)",
     main = "Volcano Plot")
abline(h = -log10(0.05), col = "brown", lty = 2) # Incluyo una línea de referencia además del código de color

# Creo un heatmap

significant_matrix <- assay(SE)[rownames(assay(SE)) %in% significant_peptides$ID, ]

group_order <- colData$Group[match(colnames(significant_matrix), colnames(assay(SE)))]

heatmap(significant_matrix, 
        scale = "row", 
        col = colorRampPalette(c("green", "white", "red"))(100), 
        main = "Heatmap",
        Rowv = NA,
        Colv = NA,
        labCol = group_order) # Indico que muestre el nombre de grupo al que pertecene la muestra