invisible(lapply(list(
  "stringr",
  "ggplot2",
  "Seurat"
), FUN = function(x) {
  suppressPackageStartupMessages(library(x, character.only = T))
}))

set.seed(snakemake@config[["seed"]])

sobj <- readRDS(file = snakemake@input[["seurat_object"]])

name <- "Cell.Cycle"

features <- list(S.Score = s.features, G2M.Score = g2m.features)
if (is.null(x = ctrl)) {
  ctrl <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
}
object.cc <- AddModuleScore(
  object = object, features = features,
  name = name, ctrl = ctrl
)
cc.columns <- grep(
  pattern = name, x = colnames(x = object.cc[[]]),
  value = TRUE
)
cc.scores <- object.cc[[cc.columns]]
rm(object.cc)
CheckGC()
assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores,
                                                               first = "S", second = "G2M", null = "G1") {
  if (all(scores < 0)) {
    return(null)
  } else {
    if (length(which(x = scores == max(scores))) > 1) {
      return("Undecided")
    } else {
      return(c(first, second)[which(x = scores == max(scores))])
    }
  }
})
cc.scores <- merge(
  x = cc.scores, y = data.frame(assignments),
  by = 0
)
colnames(x = cc.scores) <- c(
  "rownames", "S.Score", "G2M.Score",
  "Phase"
)
rownames(x = cc.scores) <- cc.scores$rownames
cc.scores <- cc.scores[, c("S.Score", "G2M.Score", "Phase")]
object[[colnames(x = cc.scores)]] <- cc.scores
if (set.ident) {
  object[["old.ident"]] <- Idents(object = object)
  Idents(object = object) <- "Phase"
}
return(object)




sobj = readRDS("../processed/seurat_processing_all/seurat_object.rds")
stress_sig <- read.table("../external/stress_signature_from_DB.txt",
  col.names = "gene"
)

tmp <- AddModuleScore(sobj,
  features = list(stress_sig$gene),
  name = "stress_signature"
)
summary(tmp[[]]["stress_signature1"])
library(ggplot2)
ggsave(
  plot = FeaturePlot(tmp, features = "stress_signature1"),
  filename = "stress_sig_umap.svg", device = "svg"
)

scores = tmp[["stress_signature1"]]
assignments <- apply(X = scores, MARGIN = 1, FUN = function(scores,
		first = "stressed", null = "not_stressed") {
		if (all(scores < 0)) {
				return(null)
		}
		else {
				if (length(which(x = scores == max(scores))) > 1) {
						return("Undecided")
				}
				else {
						return(first)
				}
		}
})

tmp[["is_stressed"]] <- assignments

ggsave(plot=DimPlot(tmp, group.by="is_stressed"), filename="umap_is_stressed.svg", device="svg")
table(tmp[[c("is_stressed", "sample")]])

