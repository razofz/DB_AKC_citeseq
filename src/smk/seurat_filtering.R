invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sobj <- readRDS(
  file = snakemake@input[["seurat_object"]]
)

print(length(Cells(sobj)))

sobj <- subset(sobj,
  subset = nFeature_RNA >
    snakemake@config[["filtering_params"]][["min_RNA_features"]] &
    nFeature_RNA <
      snakemake@config[["filtering_params"]][["max_RNA_features"]] &
    nCount_RNA >
      snakemake@config[["filtering_params"]][["min_RNA_counts"]] &
    nCount_RNA <
      snakemake@config[["filtering_params"]][["max_RNA_counts"]] &
    nCount_ADT >
      snakemake@config[["filtering_params"]][["min_ADT_counts"]] &
    nCount_ADT <
      snakemake@config[["filtering_params"]][["max_ADT_counts"]] &
    nCount_HTO >
      snakemake@config[["filtering_params"]][["min_HTO_counts"]] &
    nCount_HTO <
      snakemake@config[["filtering_params"]][["max_HTO_counts"]] &
    percent.mt <
      snakemake@config[["filtering_params"]][["max_mitochondrial_percentage"]]
)

print(length(Cells(sobj)))

classifications <- table(sobj$consensuscall, useNA = "ifany")
head(classifications)
classifications <- as.data.frame(classifications)
head(classifications)

classifications$colour <- c(NA, NA, 1, 2, 3, 4)
head(classifications)
classifications$colour <- as.factor(classifications$colour)

Idents(sobj) <- "consensuscall"
sobj <- subset(sobj, idents = str_c("HTO", 1:4))

saveRDS(sobj, snakemake@output[["seurat_object"]])
