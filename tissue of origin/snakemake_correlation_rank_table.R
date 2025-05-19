samples <- snakemake@input[["samples"]]

##
#print(samples)
refSample <- snakemake@params[["refSample"]]
samples_table <- read.csv(snakemake@input[["samples_file"]], sep="\t", header=TRUE)
# snakemake param für cell-line (für allFreq)-> möglicherweise für ALLE cell-lines?
##################################

proteinAtlas <- read.table(snakemake@input[["proteinAtlas"]],header=T,as.is=T,sep="\t")
rownames(proteinAtlas) <- proteinAtlas$GeneID

ndata <- proteinAtlas[,-1]
ndata[ndata == 0] <- NA
# at least 3 tissues with non-zero values
ndata <- ndata[apply(ndata,1,FUN=function(x) { length(which(is.na(x))) < length(x)-2 } ),] 
ndata[is.na(ndata)] <- 0.04
logndata <- log2(ndata)

##################################

tLabels <- read.table(snakemake@input[["labels"]],header=T,as.is=T,sep="\t",quote="\"")

fftColumns <- 29:52 # 160-222
#selFreq <- c("193","196","199")

library(gplots)

##################################
# Replace BH01 with sample you are using as reference in the rank comparison
refSample_name <- strsplit(tail(strsplit(refSample, "/")[[1]], 1), "-")[[1]][2]
#print(refSample)


filtered_samples_paths <- samples[sapply(samples, function(sample) {
  # Extraer el sample_name del path
  sample_name <- strsplit(tail(strsplit(sample, "/")[[1]], 1), "-")[[1]][2]
  
  # Verificar si el sample_name está en la columna 'sample' y si su ref_samples coincide con refSample
  sample_row <- samples_table[samples_table$sample == sample_name, ]
  if (nrow(sample_row) > 0 && sample_row$ref_samples == refSample_name) {
    return(TRUE)
  } else {
    return(FALSE)
  }
})]

#print(filtered_samples_paths)

fdata <- read.table(refSample,as.is=T,sep="\t",header=T,comment.char="~")
colnames(fdata) <- sub("X","",colnames(fdata))
rownames(fdata) <- fdata[,1]
fdata <- fdata[,c(1,rev(c(2:dim(fdata)[2])))]
logndata2 <- logndata[fdata[,1],]

selFreq <- c(colnames(fdata) > 189 & colnames(fdata) < 200)

refCorrelation <- cor(rowMeans(fdata[,selFreq]),logndata2[,order(names(logndata2))],use="pairwise.complete.obs")

pdf(sprintf(snakemake@output[["aveCorRank"]]),width=10,height=15)
all_correlation_results <- list()

for (sample in filtered_samples_paths) {
  # Leer el archivo de datos para el sample actual
  fdata <- read.table(sample, as.is = TRUE, sep = "\t", header = TRUE, comment.char = "~")
  
  # Limpiar nombres de columnas y filas
  colnames(fdata) <- sub("X", "", colnames(fdata))
  rownames(fdata) <- fdata[, 1]
  fdata <- fdata[, c(1, rev(c(2:dim(fdata)[2])))]
  
  # Filtrar y log transform
  logndata2 <- logndata[fdata[, 1], ]
  
  # Extraer el nombre del sample
  sample_name <- strsplit(tail(strsplit(sample, "/")[[1]], 1), "-")[[1]][2]
  
  # Seleccionar frecuencias dentro del rango deseado
  selFreq <- colnames(fdata) > 189 & colnames(fdata) < 200
  
  # Calcular correlación
  res <- cor(rowMeans(fdata[, selFreq]), logndata2[, order(names(logndata2))], use = "pairwise.complete.obs")
  
  # Filtrar correlaciones usando los labels de tejidos
  match <- intersect(colnames(res), gsub(' ', '', tLabels$Name))
  res_filtered <- res[, match, drop = FALSE]
  refCorrelation_filtered <- refCorrelation[, match, drop = FALSE]
  
  # Crear un data frame con los resultados
  correlation_df <- data.frame(
    sample_id = sample_name,  # Añadir sample ID
    category = tLabels$Category,
    description = tLabels$Type,
    tissue = match,
    correlation = as.numeric(res_filtered),
    rankDiff = rank(refCorrelation_filtered) - rank(res_filtered)
  )
  
  # Almacenar los resultados en la lista
  all_correlation_results[[sample_name]] <- correlation_df
  
  # Graficar los resultados
  par(mfrow = c(2, 1))
  textplot(head(correlation_df[order(correlation_df$rankDiff, decreasing = TRUE), ], 15))
  title(sprintf("By correlation rank difference: %s (vs. %s)", sample_name, refSample_name))
  textplot(head(correlation_df[order(correlation_df$correlation), ], 15))
  title(sprintf("By correlation: %s", sample_name))
}

dev.off()

# Combinar todos los resultados en un solo data frame
final_results <- do.call(rbind, all_correlation_results)

# Guardar los resultados en un archivo TSV
write.table(final_results, file = snakemake@output[["correlation_results_rank"]], sep = "\t", row.names = FALSE, quote = FALSE)

