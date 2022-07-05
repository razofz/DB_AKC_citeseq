invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])


find_markers <- function(idents, ident_1, ident_2, output_name) {
  # print(paste0("> idents: ", idents))
  # print(paste0("> ident_1: ", ident_1))
  # print(paste0("> ident_2: ", ident_2))
  # print(paste0("> output_name: ", output_name))
  # print(paste0("> snakemake@output[[output_name]]: ", snakemake@output[[output_name]]))
  print(paste0("> Finding differentially expressed genes for ", output_name, ".."))

  Idents(sobj) <- idents
  markers_df <- FindMarkers(sobj,
    test.use = "wilcox",
    ident.1 = ident_1,
    ident.2 = ident_2
  )
  markers_df <- markers_df[order(-markers_df$avg_log2FC),]
  # print(paste0("> dim(markers_df): ", dim(markers_df)))
  write.table(
    x = markers_df, file = snakemake@output[[output_name]],
    quote = F, sep = "\t"
  )
}


sobj <- readRDS(file = snakemake@input[["seurat_object"]])

find_markers(
  idents = "seurat_clusters", ident_1 = c(6, 9),
  ident_2 = c(1, 2, 4, 10), output_name = "deg_cl69_vs_rest_cl_ice"
)

find_markers(
  idents = "seurat_clusters", ident_1 = c(7), ident_2 = c(3, 5, 13),
  output_name = "deg_cl7_vs_rest_cl_37c_t"
)

find_markers(
  idents = "seurat_clusters", ident_1 = c(8), ident_2 = c(0, 11),
  output_name = "deg_cl8_vs_rest_cl_37c_no_t"
)

find_markers(
  idents = "hto", ident_1 = c("HTO2"), ident_2 = c("HTO4"),
  output_name = "deg_hto2_vs_hto4"
)

find_markers(
  idents = "hto", ident_1 = c("HTO1", "HTO3"), ident_2 = c("HTO4"),
  output_name = "deg_hto13_vs_hto4"
)

find_markers(
  idents = "hto", ident_1 = c("HTO1", "HTO3"), ident_2 = c("HTO2"),
  output_name = "deg_hto13_vs_hto2"
)

find_markers(
  idents = "incubation_method", ident_1 = "ice", ident_2 = "37c",
  output_name = "deg_ice_vs_37c"
)

