process_dataset_and_get_mds <- function(counts_df, groups, data_type) {
  if ("transcript_id" %in% colnames(counts_df)) {
    counts_data <- counts_df %>% dplyr::select(-transcript_id)
  } else {
    counts_data <- counts_df
  }
  if ("gene_id" %in% colnames(counts_data)) {
    counts_data <- counts_data %>% dplyr::select(-gene_id)
  }

  dgelist <- DGEList(counts = counts_data, group = groups)
  grp <- sapply(strsplit(colnames(dgelist), "-", fixed = TRUE), .subset, 1)
  mm <- model.matrix(~ 0 + grp)
  dgelist <- calcNormFactors(dgelist)
  dgelist <- estimateDisp(dgelist, mm)
  tmp <- plotMDS(dgelist, plot = FALSE)

  data.frame(
    Sample = colnames(counts_data),
    MDS1 = tmp$x,
    MDS2 = tmp$y,
    type = data_type
  )
}

plot_figure_01_d <- function(output_pdf, output_svg) {
  kinnex_tx <- vroom::vroom("results/quantify_bulk_downsampled/format/run_isosceles_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv")
  ont_tx <- vroom::vroom("results/quantify_bulk_downsampled/format/run_isosceles_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv")
  illumina_tx <- vroom::vroom("results/quantify_bulk_downsampled/format/run_salmon_illumina_corrected/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv")
  kinnex_pb_tx <- vroom::vroom("results/format_quantify_subsampled/gencode/pb/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv")
  ont_pb_tx <- vroom::vroom("results/format_quantify_subsampled/gencode/ont/1_60000000.0/oarfish/pseudobulk/transcript_counts_formatted.tsv")

  sample_groups_new <- c(0, 0, 0, 1, 1, 1)
  set.seed(42)

  df_kinnex_bulk <- process_dataset_and_get_mds(kinnex_tx, sample_groups_new, "PB (Bulk)")
  df_ont_bulk <- process_dataset_and_get_mds(ont_tx, sample_groups_new, "ONT (Bulk)")
  df_illumina_bulk <- process_dataset_and_get_mds(illumina_tx, sample_groups_new, "Illmn (Bulk)")
  df_kinnex_pb <- process_dataset_and_get_mds(kinnex_pb_tx, sample_groups_new, "PB (SC)")
  df_ont_pb <- process_dataset_and_get_mds(ont_pb_tx, sample_groups_new, "ONT (SC)")

  df_overall <- rbind(df_kinnex_bulk, df_ont_bulk, df_illumina_bulk, df_kinnex_pb, df_ont_pb)
  df_overall$Sample <- c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")
  df_overall$type <- factor(
    df_overall$type,
    levels = c("Illmn (Bulk)", "ONT (Bulk)", "ONT (SC)", "PB (Bulk)", "PB (SC)")
  )

  panel_d <- df_overall %>%
    ggplot(aes(x = MDS1, y = MDS2, label = Sample)) +
    geom_point(size = 4) +
    geom_label_repel(
      aes(label = Sample),
      min.segment.length = 0,
      size = 3,
      max.overlaps = 15,
      show.legend = FALSE
    ) +
    labs(x = "MDS1", y = "MDS2", color = "Condition", shape = "Data type") +
    facet_wrap(~ type, nrow = 1) +
    theme_big_simple()

  ggsave(output_pdf, panel_d, dpi = 600, width = 18, height = 6)
  ggsave(output_svg, panel_d, dpi = 600, width = 18, height = 6)
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
  library(edgeR)
  library(ggrepel)
  library(ggpubfigs)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_01/panel_D.pdf"
  output_svg <- "results/figure_01/panel_D.svg"
}

status <- plot_figure_01_d(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}
