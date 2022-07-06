invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sobj <- readRDS(
  file = snakemake@output[["seurat_object_demultiplexed"]]
)

ggplot(classifications, aes(x = Var1, y = Freq, fill = colour)) +
  geom_bar(stat = "identity") + scale_fill_discrete() +
  NoLegend() + xlab("Global classification") +
  ylab("barcodes")

ggsave(snakemake@output[["hto_distribution_cellhashr"]], device = "svg")
