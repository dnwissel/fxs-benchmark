calculate_dge <- function(counts,
                          gene_name,
                          filter_values = NULL,
                          cpm_threshold = 1,
                          cpm_sample_threshold = 3) {
  y <- DGEList(
    counts = counts[, -1],
    genes = data.frame(gene_id = counts$gene_id, gene_name = gene_name)
  )
  y <- calcNormFactors(y)
  if (!is.null(filter_values)) {
    filter_mat <- filter_values %>%
      dplyr::select(gene_id, all_of(colnames(y$counts))) %>%
      tibble::column_to_rownames("gene_id")
    filter_mat <- filter_mat[y$genes$gene_id, , drop = FALSE]
    keep <- rowSums(filter_mat > cpm_threshold, na.rm = TRUE) >= cpm_sample_threshold
  } else {
    keep <- rowSums(edgeR::cpm(y) > cpm_threshold) >= cpm_sample_threshold
  }

  y <- y[keep, ]
  fxs <- c(rep(0, 3), rep(1, 3))
  mm <- model.matrix(~ fxs)
  dge <- edgeR::estimateDisp(y, mm)
  fit <- edgeR::glmQLFit(dge, mm)
  mc <- limma::makeContrasts(fxs, levels = colnames(fit$coefficients))
  qlf <- edgeR::glmQLFTest(fit, contrast = mc)
  tt <- edgeR::topTags(qlf, n = Inf)
  tt$table$gene_name[is.na(tt$table$gene_name)] <- ""

  plt_frame <- data.frame(
    coef = tt$table$logFC,
    p_val = tt$table$FDR,
    gene_id = tt$table$gene_id,
    gene_name = tt$table$gene_name
  )
  plt_frame$gene_type <- ifelse(
    plt_frame$p_val < 0.1 & plt_frame$coef > 1,
    "Up",
    ifelse(plt_frame$p_val < 0.1 & plt_frame$coef < -1, "Down", "NS")
  )
  plt_frame
}

plot_figure_01_c <- function(output_pdf, output_svg) {
  kinnex <- vroom::vroom("results/quantify_bulk_downsampled/format/run_isosceles_lr/1_15000000.0/pb/gencode/gene_counts_formatted.tsv")
  ont <- vroom::vroom("results/quantify_bulk_downsampled/format/run_isosceles_lr/1_15000000.0/ont/gencode/gene_counts_formatted.tsv")
  illumina <- vroom::vroom("results/quantify_bulk_downsampled/format/run_salmon_illumina_corrected/1_15000000.0/illumina/gencode/gene_counts_formatted.tsv")
  illumina_tpm <- vroom::vroom("results/quantify_bulk_downsampled/format/run_salmon_illumina_tpm/1_15000000.0/illumina/gencode/gene_counts_formatted.tsv")
  kinnex_pseudobulk <- vroom::vroom("results/format_quantify_subsampled/gencode/pb/1_60000000.0/oarfish/pseudobulk/gene_counts_formatted.tsv")
  ont_pseudobulk <- vroom::vroom("results/format_quantify_subsampled/gencode/ont/1_60000000.0/oarfish/pseudobulk/gene_counts_formatted.tsv")

  mart <- useEnsembl("ensembl", "hsapiens_gene_ensembl")
  z <- getBM(
    c("ensembl_gene_id", "hgnc_symbol"),
    "ensembl_gene_id",
    sapply(strsplit(kinnex$gene_id, "\\."), function(x) x[[1]]),
    mart
  )
  z$hgnc_symbol[z$hgnc_symbol == ""] <- z$ensembl_gene_id[z$hgnc_symbol == ""]

  gene_name <- data.frame(gene_id = sapply(strsplit(kinnex$gene_id, "\\."), function(x) x[[1]])) %>%
    left_join(data.frame(gene_id = z$ensembl_gene_id, gene_name = z$hgnc_symbol), multiple = "first") %>%
    pull(gene_name)

  plt_frame <- rbind(
    cbind(calculate_dge(illumina, gene_name, filter_values = illumina_tpm), level = "Illmn (Bulk)"),
    cbind(calculate_dge(ont, gene_name), level = "ONT (Bulk)"),
    cbind(calculate_dge(kinnex, gene_name), level = "PB (Bulk)"),
    cbind(calculate_dge(ont_pseudobulk, gene_name), level = "ONT (SC)"),
    cbind(calculate_dge(kinnex_pseudobulk, gene_name), level = "PB (SC)")
  )

  sig_il_genes <- plt_frame %>% filter(gene_name %in% c("FMR1"))
  cols <- c("Up" = "#CC79A7", "Down" = "#56B4E9", "NS" = "grey")

  panel_c <- plt_frame %>%
    ggplot(aes(x = coef, y = -log10(p_val))) +
    geom_point(aes(colour = gene_type), alpha = 0.7, shape = 16, size = 4) +
    geom_hline(yintercept = -log10(0.1), linetype = "dashed") +
    geom_vline(xintercept = c(log2(0.5), log2(2)), linetype = "dashed") +
    facet_wrap(~ level, nrow = 1) +
    xlim(-16, 16) +
    geom_label_repel(
      data = sig_il_genes,
      aes(label = gene_name, fill = NULL),
      force = 2,
      size = 7,
      nudge_y = 1,
      segment.colour = "grey50",
      segment.size = 0.5
    ) +
    scale_color_manual(values = cols) +
    ggpubfigs::theme_big_simple() +
    labs(x = "logFC (isoB11 vs E3)", y = "-log10(p)", color = "")

  ggsave(output_pdf, panel_c, dpi = 600, width = 18, height = 6)
  ggsave(output_svg, panel_c, dpi = 600, width = 18, height = 6)
  return(0)
}

if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(vroom)
  library(biomaRt)
  library(edgeR)
  library(limma)
  library(ggrepel)
  library(ggpubfigs)
  library(tibble)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_01/panel_C.pdf"
  output_svg <- "results/figure_01/panel_C.svg"
}

status <- plot_figure_01_c(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}
