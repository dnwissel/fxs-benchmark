log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(vroom)
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
  library(scuttle)
  library(tibble)
})

t2g_path <- snakemake@input[["t2g"]]
method   <- snakemake@wildcards[["method"]]
mtx_files <- snakemake@input[["mtx"]]

process_mtx <- function(mtx_path, sample_id, method, t2g_path) {
  dir <- dirname(mtx_path)
  
  if (tolower(method) == "kallisto") {
    bc_path <- file.path(dir, "count.barcodes.txt")
    tx_path <- mtx_path
    ft_path <- file.path(dir, "transcripts.txt")
  } else if (tolower(method) == "oarfish") {
    bc_path <- gsub(".count.mtx", ".barcodes.txt", mtx_path)
    tx_path <- mtx_path
    ft_path <- gsub(".count.mtx", ".features.txt", mtx_path)
  } else {
    bc_path <- file.path(dir, "barcodes.txt")
    tx_path <- mtx_path
    ft_path <- file.path(dir, "features.txt")
  }

  t2g_map <- read.table(t2g_path, header = FALSE)
  colnames(t2g_map) <- c("transcript_id", "gene_id")
  
  tx_sparse_count_matrix <- Matrix::readMM(tx_path)
  features <- vroom::vroom(ft_path, delim = "\t", col_names = FALSE, show_col_types = FALSE)
  barcodes <- vroom::vroom(bc_path, delim = "\t", col_names = FALSE, show_col_types = FALSE)

  if (nrow(tx_sparse_count_matrix) != nrow(features)) {
    tx_sparse_count_matrix <- t(tx_sparse_count_matrix)
  }

  tx_sce <- SingleCellExperiment(
    assays = list(counts = tx_sparse_count_matrix),
    colData = data.frame(sample_id = rep(sample_id, ncol(tx_sparse_count_matrix))),
    rowData = DataFrame(t2g_map[match(features$X1, t2g_map$transcript_id), ])
  )
  
  colnames(tx_sce) <- barcodes$X1
  rownames(tx_sce) <- features$X1

  gene_sce <- aggregateAcrossFeatures(tx_sce, rowData(tx_sce)[["gene_id"]])
  assays(gene_sce)[["counts"]] <- as(assays(gene_sce)[["counts"]], "dgCMatrix")

  tx_pb <- setNames(data.frame(rowSums(counts(tx_sce)), row.names = rownames(tx_sce)), sample_id)
  gene_pb <- setNames(data.frame(rowSums(counts(gene_sce)), row.names = rownames(gene_sce)), sample_id)

  return(list(
    tx_sce = tx_sce, 
    gene_sce = gene_sce,
    tx_pb = tx_pb,
    gene_pb = gene_pb
  ))
}

results <- lapply(mtx_files, function(f) {
  if (tolower(method) == "oarfish") {
    sample_id <- gsub(".count.mtx", "", basename(f))
  } else {
    sample_id <- basename(dirname(f))
  }
  
  process_mtx(f, sample_id, method, t2g_path)
})

results <- results[!sapply(results, is.null)]


saveRDS(Reduce(cbind, lapply(results, `[[`, "tx_sce")), snakemake@output[["tx_sce"]])
saveRDS(Reduce(cbind, lapply(results, `[[`, "gene_sce")), snakemake@output[["gene_sce"]])

tx_pb_merged <- Reduce(cbind, lapply(results, `[[`, "tx_pb"))
gene_pb_merged <- Reduce(cbind, lapply(results, `[[`, "gene_pb"))

tx_pb_merged <- tx_pb_merged %>% rownames_to_column("transcript_id")
gene_pb_merged <- gene_pb_merged %>% rownames_to_column("gene_id")

write.table(tx_pb_merged, file = snakemake@output[["tx_pb"]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(gene_pb_merged, file = snakemake@output[["gene_pb"]], sep = "\t", quote = FALSE, row.names = FALSE)

sink()
sink()