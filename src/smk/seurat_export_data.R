invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sobj <- readRDS(file = snakemake@input[["seurat_object"]])

umap <- sobj[["umap"]]@cell.embeddings
write.table(umap, snakemake@output[["umap_embeddings"]], quote = F, sep = ",")

metadata <- sobj[[]]
metadata[["barcodes"]] <- rownames(metadata)
metadata <- metadata[, c(length(metadata), 1:(length(metadata) - 1))]
write.table(
  metadata,
  snakemake@output[["metadata"]],
  quote = F, sep = "\t",
  row.names = F
)
