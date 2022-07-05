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

for (modality in c("RNA", "ADT", "HTO")) {
  p <- VlnPlot(sobj,
    assay = "RNA",
    features = c(
      str_c("nCount_", modality),
      str_c("nFeature_", modality)
    )
  )
  ggsave(
    file = snakemake@output[[str_c(modality, "_vln")]],
    plot = p,
    device = "svg"
  )
  p <- VlnPlot(sobj,
    assay = "RNA",
    log = T,
    features = c(
      str_c("nCount_", modality),
      str_c("nFeature_", modality)
    )
  )
  ggsave(snakemake@output[[str_c(modality, "_vln_log")]],
    plot = p,
    device = "svg"
  )
}

sobj[["percent.mt"]] <- PercentageFeatureSet(sobj,
  pattern = "^mt-", assay = "RNA"
)
ggsave(snakemake@output[["mito_vln"]],
  plot = VlnPlot(sobj,
    assay = "RNA",
    features = "percent.mt"
  ),
  device = "svg"
)
ggsave(snakemake@output[["mito_vln_log"]],
  plot = VlnPlot(sobj,
    assay = "RNA",
    log = T,
    features = "percent.mt"
  ),
  device = "svg"
)
saveRDS(sobj, snakemake@output[["seurat_object"]])

