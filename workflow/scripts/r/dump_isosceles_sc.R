library(Matrix)
library(SingleCellExperiment)

sce_path    <- snakemake@input[["rds"]]
t2g_path    <- snakemake@input[["t2g"]]
output_path <- dirname(snakemake@output[["mtx"]])

dump_isosceles_sc <- function(sce_path, t2g_path, out_dir) {
  sce <- readRDS(sce_path)
  sce <- sce[!grepl(",", rowData(sce)$compatible_tx), ]
  
  cts <- assays(sce)[["counts"]]
  bcs <- gsub("sample.", "", colnames(sce))
  rownames(cts) <- rowData(sce)$compatible_tx
  colnames(cts) <- bcs
  
  t2g_map <- read.table(t2g_path)
  missing_tx <- setdiff(t2g_map$V1, rownames(cts))
  
  missing_mtx <- Matrix(0, 
                        nrow = length(missing_tx), 
                        ncol = ncol(cts), 
                        sparse = TRUE, 
                        dimnames = list(missing_tx, bcs))
  
  complete_cts <- rbind(cts, missing_mtx)
  write.table(rownames(complete_cts), file = file.path(out_dir, "features.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(bcs, file = file.path(out_dir, "barcodes.txt"), 
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  writeMM(t(complete_cts), file = file.path(out_dir, "count.mtx"))
}

dump_isosceles_sc(sce_path, t2g_path, output_path)
