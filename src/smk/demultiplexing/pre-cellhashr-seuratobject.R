lapply(list(
  "stringr",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
})


set.seed(snakemake@config[["seed"]])

h5 <- Read10X_h5(snakemake@input[["filtered_h5"]])

adts <- h5$`Antibody Capture`[1:4, ]
htos <- h5$`Antibody Capture`[5:8, ]
sobj <- CreateSeuratObject(h5$`Gene Expression`)
sobj[["ADT"]] <- CreateAssayObject(adts)
sobj[["HTO"]] <- CreateAssayObject(htos)

cells <- Cells(sobj)
cells <- unlist(lapply(cells, FUN = function(x) gsub("([ACGT]+)-1", "\\1", x)))

write.table(cells,
  file = snakemake@output[["cell_whitelist"]], sep = "\t", col.names = F,
  row.names = F, quote = F
)

saveRDS(sobj, file = snakemake@output[["seurat_object_unprocessed"]])
