library(cellhashR)


set.seed(snakemake@config[["seed"]])

pdf(snakemake@output[["rplots_pdf"]])

raw_count_dir <- dirname(snakemake@input[["mtx"]])

barcode_data <- ProcessCountMatrix(
  rawCountData = raw_count_dir,
  metricsFile = snakemake@output[["metrics_file_ProcessCountMatrix"]],
  barcodeBlacklist = c("unmapped")
)

methods <- c("bff_cluster", "gmm_demux", "multiseq", "htodemux", "dropletutils")

# df <- GenerateCellHashingCalls(
#   barcodeMatrix = barcodeData, methods = methods,
#   doTSNE = T, metricsFile = "hashcall_metrics.file"
# )

# write.table(df,
#   file = "hashing_results.txt", sep = "\t", row.names = FALSE,
#   quote = FALSE
# )

whitelist <- read.table(snakemake@input[["cell_whitelist"]])

df <- GenerateCellHashingCalls(
  barcodeMatrix = barcode_data, methods = methods,
  whitelist = whitelist, doTSNE = T,
  metricsFile = snakemake@output[["metrics_file_GenerateCellHashingCalls"]]
)
write.table(df,
  file = snakemake@output[["demultiplexing_results"]],
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

dev.off()
