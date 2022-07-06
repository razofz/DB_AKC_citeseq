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
  ggsave(
    file = snakemake@output[[str_c(modality, "_vln_log")]],
    plot = p,
    device = "svg"
  )
}

ggsave(
  file = snakemake@output[["mito_vln"]],
  plot = VlnPlot(sobj,
    assay = "RNA",
    features = "percent.mt"
  ),
  device = "svg"
)
ggsave(
  file = snakemake@output[["mito_vln_log"]],
  plot = VlnPlot(sobj,
    assay = "RNA",
    log = T,
    features = "percent.mt"
  ),
  device = "svg"
)

classifications <- table(sobj$consensuscall, useNA = "ifany")
head(classifications)
classifications <- as.data.frame(classifications)
head(classifications)

classifications$colour <- c(NA, NA, 1, 2, 3, 4)
head(classifications)
classifications$colour <- as.factor(classifications$colour)

ggplot(classifications, aes(x = Var1, y = Freq, fill = colour)) +
  geom_bar(stat = "identity") + scale_fill_discrete() +
  NoLegend() + xlab("Global classification") +
  ylab("barcodes")

ggsave(snakemake@output[["hto_distribution_cellhashr"]], device = "svg")

Idents(sobj) <- "consensuscall"
sobj <- subset(sobj, idents = str_c("HTO", 1:4))

saveRDS(sobj, snakemake@output[["seurat_object"]])
