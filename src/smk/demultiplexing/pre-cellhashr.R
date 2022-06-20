library(stringr)
library(ggplot2)
library(Seurat)
library(Matrix)

set.seed(snakemake@config[["seed"]])

rawcounts <- Read10X(snakemake@input[["raw_counts"]])

raws <- rawcounts$`Antibody Capture`
raws <- raws[5:8, ]
raws

writeMM(raws, snakemake@output[["mtx"]])
write.table(colnames(raws),
  file = snakemake@output[["barcodes"]], sep = "\t",
  row.names = F, quote = F, col.names = F
)

df <- rownames(raws)
df <- as.data.frame(df)

df$short <- unlist(lapply(df$df, FUN = function(x) {
  str_split(x, pattern = "-")[[1]]
}))[c(
  1, 4,
  7, 10
)]

df <- df[, c(2, 1)]
df$type <- "Antibody Capture"
write.table(df,
  file = snakemake@output[["features"]], sep = "\t", row.names = F,
  quote = F, col.names = F
)

# ownmtx <- Read10X("antibody_mtx/")

# hashes <- read.table("cellhashr/hashing_results.txt")

# called <- hashes[str_starts(hashes$consensuscall, "HTO"), ]

# rownames(called) <- str_c(rownames(called), "-1")
# length(intersect(Cells(sobj), rownames(called)))
# rownames(hashes) = str_c(rownames(hashes), "-1")
# length(intersect(Cells(sobj), rownames(hashes)))
# length(difference(Cells(sobj), rownames(hashes)))
# length(setdiff(Cells(sobj), rownames(hashes)))
# length(setdiff(rownames(hashes), Cells(sobj)))

