invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sobj <- readRDS(file = snakemake@input[["seurat_object"]])

sobj[["orig.ident"]] <- snakemake@config[["project_name"]]
sobj[["consensuscall.global"]] <- NULL
sobj[["hto"]] <- sobj[["consensuscall"]]
sobj[["consensuscall"]] <- NULL
sobj[["sample"]] <- sobj[["hto"]]
sobj[["sample_longname"]] <- sobj[["hto"]]
sobj[["buffer_treatment"]] <- sobj[["hto"]]
sobj[["incubation_method"]] <- sobj[["hto"]]

for (i in seq_len(length(snakemake@config[["samples"]]))) {
  sobj[["sample"]][sobj[["sample"]] == str_c("HTO", i)] <-
    snakemake@config[["samples"]][[i]]
  sobj[["sample_longname"]][sobj[["sample_longname"]] == str_c("HTO", i)] <-
    snakemake@config[["samples_longname"]][[i]]
  if (i %in% snakemake@config[["sample_groups"]][["with_triptolide"]]) {
    sobj[["buffer_treatment"]][sobj[["buffer_treatment"]] == str_c("HTO", i)] <-
      "with_triptolide"
  } else if (i %in% snakemake@config[["sample_groups"]][["no_triptolide"]]) {
    sobj[["buffer_treatment"]][sobj[["buffer_treatment"]] == str_c("HTO", i)] <-
      "no_triptolide"
  }
  if (i %in% snakemake@config[["sample_groups"]][["37c"]]) {
    sobj[["incubation_method"]][
      sobj[["incubation_method"]] == str_c("HTO", i)
    ] <-
      "37c"
  } else if (i %in% snakemake@config[["sample_groups"]][["ice"]]) {
    sobj[["incubation_method"]][
      sobj[["incubation_method"]] == str_c("HTO", i)
    ] <-
      "ice"
  }
}

sobj <- NormalizeData(sobj, assay = "RNA")
sobj <- FindVariableFeatures(sobj, assay = "RNA")
write.table(VariableFeatures(sobj), snakemake@output[["hvgs"]],
  quote = F, row.names = F, col.names = F
)
sobj <- ScaleData(sobj,
  assay = "RNA",
  features = VariableFeatures(sobj, assay = "RNA")
)
sobj <- NormalizeData(sobj, assay = "HTO")
sobj <- NormalizeData(sobj, assay = "ADT")

sobj <- CellCycleScoring(sobj,
  s.features = str_to_title(cc.genes$s.genes),
  g2m.features = str_to_title(cc.genes$g2m.genes),
  set.ident = T
)
# sobj <- ScaleData(sobj,
#   vars.to.regress = c("S.Score", "G2M.Score"),
#   assay = "RNA",
#   features = rownames(sobj)
# )

sobj <- RunPCA(sobj, assay = "RNA", seed.use = snakemake@config[["seed"]])
sobj <- FindNeighbors(sobj, assay = "RNA", reduction = "pca", dims = 1:30)
sobj <- FindClusters(sobj, graph.name = "RNA_snn", resolution = .8)
sobj <- RunUMAP(sobj,
  dims = 1:30, reduction = "pca",
  assay = "RNA", return.model = T
)

saveRDS(sobj, snakemake@output[["seurat_object"]])
