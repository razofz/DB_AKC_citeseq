invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sobj <- readRDS(file = snakemake@input[["seurat_object_unprocessed"]])

hasheswhitelist <- read.table(snakemake@input[["demultiplexing_results"]])
rownames(hasheswhitelist) <- str_c(rownames(hasheswhitelist), "-1")

length(setdiff(Cells(sobj), rownames(hasheswhitelist)))

subsetted <- hasheswhitelist[Cells(sobj), c(
  "consensuscall",
  "consensuscall.global"
)]

table(subsetted$consensuscall, useNA = "ifany")
table(rownames(subsetted) == rownames(sobj[[]]))

sobj[["consensuscall"]] <- subsetted[["consensuscall"]]
sobj[["consensuscall.global"]] <- subsetted[["consensuscall.global"]]

classifications <- table(sobj$consensuscall, useNA = "ifany")
classifications <- as.data.frame(classifications)

classifications$colour <- c(NA, NA, 1, 2, 3, 4, NA, NA)
classifications$colour <- as.factor(classifications$colour)

sobj[["percent.mt"]] <- PercentageFeatureSet(sobj,
  pattern = "^mt-", assay = "RNA"
)

saveRDS(sobj, snakemake@output[["seurat_object_demultiplexed"]])
