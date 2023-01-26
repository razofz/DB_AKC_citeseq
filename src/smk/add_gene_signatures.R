invisible(lapply(list(
  "stringr",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

################################################################################
#                                  Load data                                   #
################################################################################

sobj <- readRDS(snakemake@input[["seurat_object"]])

DefaultAssay(sobj) <- "RNA"

stress_signature <- read.delim(snakemake@input[["stress_signature"]])

################################################################################
#                       Function to add signature score                        #
################################################################################

# Slightly modified version of the Seurat::CellCycleScoring method
add_signature_score <- function(seurat_object, signature, signature_name) {
  sig_col <- str_c(signature_name, "_signature")
  signature_name <- str_c(signature_name, "ed")
  result <- AddModuleScore(seurat_object,
    features = signature, name = sig_col
  )
  summary(result[[]][str_c(sig_col, "1")])
  scores <- result[[str_c(sig_col, "1")]]
  assignments <- apply(
    X = scores, MARGIN = 1,
    FUN = function(scores, first = signature_name, null = stringr::str_c(
                     "not_",
                     signature_name
                   )) {
      if (all(scores < 0)) {
        return(null)
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return("Undecided")
        } else {
          return(first)
        }
      }
    }
  )
  result[[str_c("is_", signature_name)]] <- assignments
  return(result)
}

################################################################################
#          Use the function to add scores for the stress signature             #
################################################################################

sobj <- add_signature_score(sobj, stress_signature, "stress")

################################################################################
#                            Save the seurat object                            #
################################################################################

saveRDS(sobj, snakemake@output[["seurat_object"]])

# DimPlot(sobj, group.by = "is_HSC") + coord_fixed()
# DimPlot(sobj, group.by = "is_earlydiff") + coord_fixed()
# table(sobj[[c("seurat_clusters", "is_HSC")]])
# df <- sobj[[c("seurat_clusters", "is_HSC", "is_earlydiff")]]
# table(df[df$seurat_cluster == "0", ][c("is_HSC", "is_earlydiff")])

# VlnPlot(sobj, features = c("nCount_RNA", "nFeature_RNA"), log = T)
# FeaturePlot(sobj, features = "Cd48") + coord_fixed()
# cells_hsc_nodiff <- WhichCells(sobj, expression = ((is_HSC == "HSC") &
#   (is_earlydiff == "not_earlydiff")))
# length(cells_hsc_nodiff)

# DimPlot(sobj, cells.highlight = cells_hsc_nodiff) + coord_fixed()

# sobj[[]][c(sobj[[]]$is_HSC == "HSC" &
#   sobj[[]]$is_earlydiff == "not_earlydiff"), ]
