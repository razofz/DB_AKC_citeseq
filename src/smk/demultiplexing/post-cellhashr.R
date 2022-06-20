library(stringr)
library(ggplot2)
library(Seurat)


set.seed(snakemake@config[["seed"]])

sobj <- readRDS(file = snakemake@input[["seurat_object_unprocessed"]])

hasheswhitelist <- read.table(snakemake@input[["demultiplexing_results"]])
rownames(hasheswhitelist) <- str_c(rownames(hasheswhitelist), "-1")

length(setdiff(Cells(sobj), rownames(hasheswhitelist)))

sobj[["consensuscall"]] <- list()

subsetted <- hasheswhitelist[Cells(sobj), c(
  "consensuscall",
  "consensuscall.global"
)]

table(subsetted$consensuscall, useNA = "ifany")
table(rownames(subsetted) == rownames(sobj[[]]))

sobj[["consensuscall"]] <- subsetted[["consensuscall.global"]]
sobj[["consensuscall.global"]] <- subsetted[["consensuscall.global"]]

table(sobj[[c("HTO_classification.global", "consensuscall.global")]])

classifications <- table(sobj$consensuscall, useNA = "ifany")
classifications <- as.data.frame(classifications)

classifications$colour <- c(NA, NA, 1, 2, 3, 4, NA, NA)
classifications$colour <- as.integer(classifications$colour)

ggplot(classifications, aes(x = Var1, y = Freq, fill = colour)) +
  geom_bar(stat = "identity") + scale_fill_discrete() +
  NoLegend() + xlab("Global classification") +
  ylab("barcodes")

ggsave(snakemake@output[["hto_distribution_cellhashr"]], device = "svg")
saveRDS(sobj, snakemake@output[["seurat_object_demultiplexed"]])
