invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sobj <- readRDS(file = snakemake@input[["seurat_object"]])

p <- DimPlot(sobj,
  group.by = "seurat_clusters", label = T, label.size = 2,
  label.box = T) + coord_fixed()
ggsave(snakemake@output[["umap_clusters"]], plot = p, device = "svg")
p <- DimPlot(sobj, group.by = "seurat_clusters") + coord_fixed()
ggsave(snakemake@output[["umap_clusters_unlabelled"]], plot = p, device = "svg")
p <- DimPlot(sobj,
  group.by = "sample", label = T, label.size = 3, repel = T,
  label.box = T
) + coord_fixed()
ggsave(snakemake@output[["umap_sample"]], plot = p, device = "svg")
p <- DimPlot(sobj,
  group.by = "hto", label = T, label.size = 3, repel = T,
  label.box = T
) + coord_fixed()
ggsave(snakemake@output[["umap_hto"]], plot = p, device = "svg")
p <- DimPlot(sobj,
  group.by = "incubation_method", label = T, repel = T,
  label.box = T
) + coord_fixed()
ggsave(snakemake@output[["umap_incubation_method"]], plot = p, device = "svg")
p <- DimPlot(sobj, group.by = "buffer_treatment") + coord_fixed()
ggsave(snakemake@output[["umap_buffer_treatment"]], plot = p, device = "svg")
p <- DimPlot(sobj, group.by = "Phase") + coord_fixed()
ggsave(snakemake@output[["umap_cellcycle"]], plot = p, device = "svg")

