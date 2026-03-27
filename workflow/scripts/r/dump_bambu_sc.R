library(Matrix)
library(SingleCellExperiment)

sce_path    <- snakemake@input[["rds"]]
sample_id   <- snakemake@wildcards[["sample"]]
output_path <- dirname(snakemake@output[["mtx"]])

dump_bambu_sc <- function(sce_path, id, out_dir) {
  sce <- readRDS(sce_path)
  cts <- assays(sce)[["counts"]]
  
  bcs <- gsub(paste0(id, "_"), "", colnames(sce))
  colnames(cts) <- bcs
  
  write.table(rownames(cts), file = file.path(out_dir, "features.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(bcs, file = file.path(out_dir, "barcodes.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeMM(t(cts), file = file.path(out_dir, "count.mtx"))
}

dump_bambu_sc(sce_path, sample_id, output_path)
