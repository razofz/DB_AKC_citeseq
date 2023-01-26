invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sobj <- readRDS(file = snakemake@input[["seurat_object"]])

p <- DimPlot(sobj, group.by = "sample", label = F) + coord_fixed()
ggsave(snakemake@output[["umap_sample"]], plot = p, device = "svg")

p <- DimPlot(sobj, group.by = "is_stressed", label = F) + coord_fixed()
ggsave(snakemake@output[["umap_stress_sig_classification"]],
  plot = p, device = "svg"
)

# p <- FeaturePlot(sobj, features = "stress_signature1") + coord_fixed()
# ggsave(snakemake@output[["umap_stress_sig_score"]], plot = p, device = "svg")

p <- FeaturePlot(sobj,
  features = "stress_signature1", cols = c(
    "gray", "lightcoral", "red"
  )
) + coord_fixed()
ggsave(snakemake@output[["umap_stress_sig_score_q100"]],
  plot = p, device = "svg"
)

p <- FeaturePlot(sobj,
  features = "stress_signature1", cols = c(
    "gray", "lightcoral", "red"
  ), max.cutoff = "q99"
) + coord_fixed()
ggsave(snakemake@output[["umap_stress_sig_score_q99"]],
  plot = p, device = "svg"
)

p <- FeaturePlot(sobj,
  features = "stress_signature1", cols = c(
    "lightblue", "red"
  )
) + coord_fixed()
ggsave(snakemake@output[["umap_stress_sig_score_blue"]],
  plot = p, device = "svg"
)

