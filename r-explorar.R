library(SummarizedExperiment)
se_cachexia <- readRDS("ruta/al/archivo")

print("Resumen:")
summary(se_cachexia)

print("Estructura:")
str(se_cachexia)

print("Matriz de expresión:")
head(assay(se_cachexia))

print("Metadatos:")
head(rowData(se_cachexia))

print("Generando boxplot de los datos de expresión:")
boxplot(assay(se_cachexia))