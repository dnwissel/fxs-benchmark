log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(dplyr)
  library(harmony)
  library(scater)
  library(scDblFinder)
  library(scran)
  library(scuttle)
  library(SingleCellExperiment)
})


filter_sce <- function(sce_path,
                       gtf_path,
                       tot_counts_filter,
                       mito_pct_filter,
                       n_detected_genes_filter,
                       tech,
                       mgv_bio_t,
                       mgv_fdr_t,
                       mgv_mean_t,
                       sample_specific_cutoff = TRUE) {

  # Load data
  sce <- readRDS(sce_path)
  
  # Keep spliced counts
  sce <- sce[!grepl("-I$", rownames(sce)), ]
  
  # Add gene_name to rowData
  gtf <- rtracklayer::import(gtf_path)
  row_meta <- gtf %>% 
    as.data.frame() %>% 
    distinct(gene_id, gene_name)
  
  rd <- rowData(sce) %>% as.data.frame %>%
    left_join(row_meta, multiple = "first")
  rowData(sce) <- rd
  
  # Find mitochondrial genes and add QC metrics
  grep("^MT-", rd$gene_name, value = TRUE)
  sce <- addPerCellQCMetrics(sce, subsets=list(mito=grepl("^MT-", rd$gene_name)))
  
  ## Plots before filtering
  cd <- colData(sce) %>% as.data.frame
  ggplot(cd, aes(x = total,
                 y = detected,
                 colour = subsets_mito_percent)) +
    geom_point() +
    scale_x_log10() + scale_y_log10() +
    geom_density_2d(colour = "orange") +
    facet_wrap(~sample_id) +
    theme_minimal() +
    ggtitle("#detected genes vs #total counts")
  
  ggplot(cd, aes(x = total, y = subsets_mito_percent)) +
    geom_point() + scale_x_log10() + scale_y_sqrt() +
    facet_wrap(~sample_id, nrow = 2) + 
    geom_hline(yintercept=c(1, mito_pct_filter), colour="orange") +
    geom_vline(xintercept=c(tot_counts_filter), colour="orange") +
    geom_density2d() +
    theme_minimal() +
    ggtitle("#total counts vs mitochondrial genes %")
  
  ggplot(cd, aes(x = detected, y = subsets_mito_percent)) +
    geom_point() + scale_x_log10() + scale_y_sqrt() +
    facet_wrap(~sample_id, nrow = 2) + 
    geom_hline(yintercept=c(1, mito_pct_filter), colour="orange") +
    geom_vline(xintercept=c(n_detected_genes_filter), colour="orange") +
    geom_density2d() +
    theme_minimal() +
    ggtitle("#detected genes vs mitochondrial genes %")
  
  # Filter single-cell object
  sce <- sce[, !is.na(sce$subsets_mito_percent)]
  sce <- sce[,sce$subsets_mito_percent<mito_pct_filter]
  sce <- sce[,sce$total>tot_counts_filter]
  
  if(sample_specific_cutoff) {
    # Find sample-specific cutoffs for counts
    totals <- split(log10(sce$total), sce$sample_id)
    rng <- range(unlist(totals))
    
    densities <- mapply(function(u, v) {
      d <- density(u, bw = 0.05, from = rng[1], to = rng[2])
      data.frame(sample_id = v, x = d$x, y = d$y)
    }, totals, names(totals), SIMPLIFY = FALSE)
    
    cutoffs <- sapply(densities, function(u) {
      n <- nrow(u)
      dx <- diff(u$y)
      u$x[dplyr::first(which(dx[-1] > 0 & dx[-(n-1)] < 0))]+diff(u$x)[1]/2
    })
    
    par(mfrow=c(1,3))
    for(i in 1:length(densities)) {
      plot(densities[[i]]$x, densities[[i]]$y, type="l",
           xlab="", ylab="", main=names(densities)[i])
      abline(v=cutoffs[i])
    }
    
    # Filter based on individual cutoffs to remove cells with a few counts
    sce <- sce[,sce$total > 10^cutoffs[sce$sample_id]]
  }
  
  # Compute log-normalized counts
  sce <- scuttle::logNormCounts(sce)
  
  # Find highly variable genes
  mgv <- modelGeneVar(sce, block = sce$sample_id)
  idx <- mgv$bio > mgv_bio_t & mgv$FDR < mgv_fdr_t & mgv$mean > mgv_mean_t
  mgv[idx,] %>% 
    as.data.frame %>% rownames -> hvgs
  
  # Remove mitochondrial genes from hvg, if there is any
  hvgs_mt <- rownames(grep("^MT-", rowData(sce)[hvgs,"gene_name"], value = TRUE))
  hvgs <- setdiff(hvgs, hvgs_mt)
  
  # Run harmony
  sce <- fixedPCA(sce, subset.row=hvgs, rank = 50)
  sce <- RunHarmony(
    sce,
    "sample_id"
  )
  set.seed(9)
  sce <- runUMAP(sce, dimred="HARMONY")
  
  # Find doublets
  sce <- scDblFinder(sce, samples = sce$sample_id,
                     verbose = TRUE)
  plotUMAP(sce, colour_by = "scDblFinder.score")
  
  # Filter doublets
  sce <- sce[,sce$scDblFinder.class=="singlet"]
  
  return(sce)
}

gtf_path <- snakemake@input[["gtf"]]

# ILLUMINA
tot_counts_filter <- snakemake@params[["ill_tot_counts_filter"]]
mito_pct_filter <- snakemake@params[["ill_mito_pct_filter"]]
n_detected_genes_filter <- snakemake@params[["ill_n_detected_genes_filter"]]
sce_path <- snakemake@input[["ill_gene_sce"]]
mgv_bio_t <- snakemake@params[["mgv_bio_t"]]
mgv_fdr_t <- snakemake@params[["mgv_fdr_t"]]
mgv_mean_t <- snakemake@params[["ill_mgv_mean_t"]]

sce <- filter_sce(sce_path,
                  gtf_path,
                  tot_counts_filter,
                  mito_pct_filter,
                  n_detected_genes_filter,
                  tech,
                  mgv_bio_t,
                  mgv_fdr_t,
                  mgv_mean_t)

if(!dir.exists(dirname(snakemake@output[["ill_gene_sce"]]))) dir.create(dirname(snakemake@output[["ill_gene_sce"]]), recursive=TRUE)
saveRDS(sce, snakemake@output[["ill_gene_sce"]])

# PACBIO
tot_counts_filter <- snakemake@params[["pb_tot_counts_filter"]]
mito_pct_filter <- snakemake@params[["pb_mito_pct_filter"]]
n_detected_genes_filter <- snakemake@params[["pb_n_detected_genes_filter"]]
sce_path <- snakemake@input[["pb_gene_sce"]]
mgv_bio_t <- snakemake@params[["mgv_bio_t"]]
mgv_fdr_t <- snakemake@params[["mgv_fdr_t"]]
mgv_mean_t <- snakemake@params[["pb_mgv_mean_t"]]

sce <- filter_sce(sce_path,
                  gtf_path,
                  tot_counts_filter,
                  mito_pct_filter,
                  n_detected_genes_filter,
                  tech,
                  mgv_bio_t,
                  mgv_fdr_t,
                  mgv_mean_t)

if(!dir.exists(dirname(snakemake@output[["pb_gene_sce"]]))) dir.create(dirname(snakemake@output[["pb_gene_sce"]]), recursive=TRUE)
saveRDS(sce, snakemake@output[["pb_gene_sce"]])

# ONT
tot_counts_filter <- snakemake@params[["ont_tot_counts_filter"]]
mito_pct_filter <- snakemake@params[["ont_mito_pct_filter"]]
n_detected_genes_filter <- snakemake@params[["ont_n_detected_genes_filter"]]
sce_path <- snakemake@input[["ont_gene_sce"]]
mgv_bio_t <- snakemake@params[["mgv_bio_t"]]
mgv_fdr_t <- snakemake@params[["mgv_fdr_t"]]
mgv_mean_t <- snakemake@params[["ont_mgv_mean_t"]]

sce <- filter_sce(sce_path,
                  gtf_path,
                  tot_counts_filter,
                  mito_pct_filter,
                  n_detected_genes_filter,
                  tech,
                  mgv_bio_t,
                  mgv_fdr_t,
                  mgv_mean_t)

if(!dir.exists(dirname(snakemake@output[["ont_gene_sce"]]))) dir.create(dirname(snakemake@output[["ont_gene_sce"]]), recursive=TRUE)
saveRDS(sce, snakemake@output[["ont_gene_sce"]])

sink()
sink()