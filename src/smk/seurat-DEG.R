invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])


find_markers_on_cells <- function(obj, assay, cells_1, cells_2, output_name) {
  print(paste0("> Finding differentially expressed genes for ", output_name, ".."))
  markers_df <- FindMarkers(obj[[assay]],
    test.use = "wilcox",
    cells.1 = cells_1,
    cells.2 = cells_2
  )
  markers_df <- markers_df[order(-markers_df$avg_log2FC), ]
  write.table(
    x = markers_df, file = snakemake@output[[output_name]],
    quote = F, sep = "\t"
  )
}

sobj <- readRDS(file = snakemake@input[["seurat_object"]])

cells <- list()
Idents(sobj) <- "sample"
cells$ice <- WhichCells(sobj, idents = c("ice_t", "ice_no_t"))
cells$i37c <- WhichCells(sobj, idents = c("37c_t", "37c_no_t"))
cells$i37c_t <- WhichCells(sobj, idents = "37c_t")
cells$i37c_no_t <- WhichCells(sobj, idents = "37c_no_t")
cells$ice_t <- WhichCells(sobj, idents = "ice_t")
cells$ice_no_t <- WhichCells(sobj, idents = "ice_no_t")

Idents(sobj) <- "buffer_treatment"
cells$with_triptolide <- WhichCells(sobj, idents = "with_triptolide")
cells$no_triptolide <- WhichCells(sobj, idents = "no_triptolide")

Idents(sobj) <- "Phase"
cells$S <- WhichCells(sobj, idents = "S")
cells$G1 <- WhichCells(sobj, idents = "G1")
cells$G2M <- WhichCells(sobj, idents = "G2M")

lapply(cells, length)

for (phase in c("S", "G1", "G2M")) {
  find_markers_on_cells(
    obj = sobj, assay = "RNA",
    cells_1 = intersect(cells$ice, cells[[phase]]),
    cells_2 = intersect(cells$i37c_t, cells[[phase]]),
    output_name = str_c("deg_", phase, "_ice_vs_37c_t")
  )
}

for (phase in c("S", "G1", "G2M")) {
  find_markers_on_cells(
    obj = sobj, assay = "RNA",
    cells_1 = intersect(cells$ice, cells[[phase]]),
    cells_2 = intersect(cells$i37c_no_t, cells[[phase]]),
    output_name = str_c("deg_", phase, "_ice_vs_37c_no_t")
  )
}

for (phase in c("S", "G1", "G2M")) {
  find_markers_on_cells(
    obj = sobj, assay = "RNA",
    cells_1 = intersect(cells$ice, cells[[phase]]),
    cells_2 = intersect(cells$i37c_t, cells[[phase]]),
    output_name = str_c("deg_", phase, "_ice_vs_37c_t")
  )
}

for (phase in c("S", "G1", "G2M")) {
  find_markers_on_cells(
    obj = sobj, assay = "RNA",
    cells_1 = intersect(cells$ice, cells[[phase]]),
    cells_2 = intersect(cells$i37c_t, cells[[phase]]),
    output_name = str_c("deg_", phase, "_ice_vs_37c_t")
  )
}

for (phase in c("S", "G1", "G2M")) {
  find_markers_on_cells(
    obj = sobj, assay = "RNA",
    cells_1 = intersect(cells$ice_t, cells[[phase]]),
    cells_2 = intersect(cells$i37c_t, cells[[phase]]),
    output_name = str_c("deg_", phase, "_triptolides")
  )
}

for (phase in c("S", "G1", "G2M")) {
  find_markers_on_cells(
    obj = sobj, assay = "RNA",
    cells_1 = intersect(cells$ice_no_t, cells[[phase]]),
    cells_2 = intersect(cells$i37c_no_t, cells[[phase]]),
    output_name = str_c("deg_", phase, "_no_triptolides")
  )
}

for (phase in c("S", "G1", "G2M")) {
  find_markers_on_cells(
    obj = sobj, assay = "RNA",
    cells_1 = intersect(cells$i37c_t, cells[[phase]]),
    cells_2 = intersect(cells$i37c_no_t, cells[[phase]]),
    output_name = str_c("deg_", phase, "_37c_no_t_vs_37c_t")
  )
