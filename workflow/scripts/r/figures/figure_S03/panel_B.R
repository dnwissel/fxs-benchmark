plot_figure_s03_b <- function(output_pdf, output_svg) {

  config <- list(log_expr_pseudo = 0.1)
  kinnex_kallisto_counts <- vroom::vroom("results/format/run_kallisto_long_lr/1_20000.0/pb/ercc/transcript_counts_formatted.tsv")
  short_sirvs <- grep("SIRV", kinnex_kallisto_counts$transcript_id, value = TRUE)
  erccs <- kinnex_kallisto_counts$transcript_id[!kinnex_kallisto_counts$transcript_id %in% short_sirvs]
  kinnex_kallisto_counts <- vroom::vroom("results/format/run_kallisto_long_lr/1_20000.0/pb/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  ont_kallisto_counts <- vroom::vroom("results/format/run_kallisto_long_lr/1_20000.0/ont/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  kinnex_bambu_counts <- vroom::vroom("results/format/run_bambu_lr/1_20000.0/pb/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  ont_bambu_counts <- vroom::vroom("results/format/run_bambu_lr/1_20000.0/ont/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  kinnex_isoquant_counts <- vroom::vroom("results/format/run_isoquant_lr/1_20000.0/pb/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  ont_isoquant_counts <- vroom::vroom("results/format/run_isoquant_lr/1_20000.0/ont/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  kinnex_oarfish_counts <- vroom::vroom("results/format/run_oarfish_lr/1_20000.0/pb/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  ont_oarfish_counts <- vroom::vroom("results/format/run_oarfish_lr/1_20000.0/ont/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  kinnex_miniquant_counts <- vroom::vroom("results/format/run_miniquant_lr/1_20000.0/pb/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  ont_miniquant_counts <- vroom::vroom("results/format/run_miniquant_lr/1_20000.0/ont/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  kinnex_isosceles_counts <- vroom::vroom("results/format/run_isosceles_lr/1_20000.0/pb/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  ont_isosceles_counts <- vroom::vroom("results/format/run_isosceles_lr/1_20000.0/ont/ercc/transcript_counts_formatted.tsv") %>% filter(transcript_id %in% erccs) %>% arrange(desc(transcript_id))
  prepped_counts <- lapply(list(ont_oarfish_counts, kinnex_oarfish_counts, ont_isoquant_counts, kinnex_isoquant_counts, ont_bambu_counts, kinnex_bambu_counts, ont_kallisto_counts, kinnex_kallisto_counts, ont_miniquant_counts, kinnex_miniquant_counts, ont_isosceles_counts, kinnex_isosceles_counts), function(x) {
      x <- x[x$transcript_id %in% erccs, ] %>% arrange(desc(transcript_id))
      x <- DGEList(x[, -(1:2)], genes = x[, 2], group = c(rep("E3", 3), rep("isoB11", 3)))
      cpm <- log2(edgeR::cpm(x, log = FALSE) + config$log_expr_pseudo)
      return(cpm)
  })
  annot <- data.frame(read_excel("results/ercc_formatted.xlsx"), check.names = FALSE)
  kinnex_kallisto_counts <- vroom::vroom("results/format/run_kallisto_long_lr/1_20000.0/pb/ercc/transcript_counts_formatted.tsv")
  short_sirvs <- grep("SIRV[0-9][0-9][0-9]$", kinnex_kallisto_counts$transcript_id, value = TRUE)
  kinnex_kallisto_counts <- vroom::vroom(glue(glue("results/format/run_kallisto_long_lr/1_75000.0/pb/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  ont_kallisto_counts <- vroom::vroom(glue(glue("results/format/run_kallisto_long_lr/1_75000.0/ont/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  kinnex_bambu_counts <- vroom::vroom(glue(glue("results/format/run_{'bambu'}_lr/1_75000.0/pb/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  ont_bambu_counts <- vroom::vroom(glue(glue("results/format/run_{'bambu'}_lr/1_75000.0/ont/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  kinnex_isoquant_counts <- vroom::vroom(glue(glue("results/format/run_{'isoquant'}_lr/1_75000.0/pb/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  ont_isoquant_counts <- vroom::vroom(glue(glue("results/format/run_{'isoquant'}_lr/1_75000.0/ont/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  kinnex_oarfish_counts <- vroom::vroom(glue(glue("results/format/run_{'oarfish'}_lr/1_75000.0/pb/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  ont_oarfish_counts <- vroom::vroom(glue(glue("results/format/run_{'oarfish'}_lr/1_75000.0/ont/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  kinnex_miniquant_counts <- vroom::vroom(glue(glue("results/format/run_{'miniquant'}_lr/1_75000.0/pb/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  ont_miniquant_counts <- vroom::vroom(glue(glue("results/format/run_{'miniquant'}_lr/1_75000.0/ont/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  kinnex_isosceles_counts <- vroom::vroom(glue(glue("results/format/run_{'isosceles'}_lr/1_75000.0/pb/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  ont_isosceles_counts <- vroom::vroom(glue(glue("results/format/run_{'isosceles'}_lr/1_75000.0/ont/sirv/transcript_counts_formatted.tsv"))) %>% filter(transcript_id %in% short_sirvs) %>% arrange(desc(transcript_id))
  prepped_counts <- lapply(list(ont_oarfish_counts, kinnex_oarfish_counts, ont_isoquant_counts, kinnex_isoquant_counts, ont_bambu_counts, kinnex_bambu_counts, ont_kallisto_counts, kinnex_kallisto_counts, ont_miniquant_counts, kinnex_miniquant_counts, ont_isosceles_counts, kinnex_isosceles_counts), function(x) {
      print(x)
      x <- x[x$transcript_id %in% short_sirvs, ] %>% arrange(desc(transcript_id))
      x <- DGEList(x[, -(1:2)], genes = x[, 2], group = c(rep("E3", 3), rep("isoB11", 3)))
      cpm <- log2(edgeR::cpm(x, log = FALSE) + config$log_expr_pseudo)
      return(cpm)
  })
  annot <- data.frame(read_excel("results/SIRV_Set1_Norm_sequence-design-overview_20210507a.xlsx"), check.names = FALSE)
  annot_df <- data.frame(transcript_id = annot$...6[-(1:3)][1:100], mw = as.numeric(annot$...8[-(1:3)][1:100]), e0 = as.numeric(annot$...10[-(1:3)][1:100]), e1 = as.numeric(annot$...11[-(1:3)][1:100]), e2 = as.numeric(annot$...12[-(1:3)][1:100])) %>% drop_na() %>% pivot_longer(!transcript_id, names_to = "spikein", values_to = "molarity")
  gtf_helper <- import(here::here("results/sirvs_gtf.gtf"))
  gtf_helper <- data.frame(gtf_helper)
  sorted_transcript_length <- gtf_helper %>% group_by(transcript_id) %>% dplyr::summarise(length = sum(width)) %>% arrange(desc(transcript_id)) %>% pull(length)
  annotation_frame_relative <- data.frame(transcript_id = gtf_helper %>% group_by(transcript_id) %>% dplyr::summarise(length = sum(width)) %>% arrange(desc(transcript_id)) %>% pull(transcript_id), len = sorted_transcript_length, molarity_fold_change = annot_df %>% filter(spikein == "e2") %>% arrange(desc(transcript_id)) %>% pull(molarity)/annot_df %>% filter(spikein == "e1") %>% arrange(desc(transcript_id)) %>% pull(molarity))
  annotation_frame_relative <- rbind(annotation_frame_relative, annotation_frame_relative, annotation_frame_relative)
  annotation_frame_relative$sample <- c(rep(1, 69), rep(2, 69), rep(3, 69))
  annotation_frame_relative <- rbind(annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative, annotation_frame_relative)
  annotation_frame_relative$tech <- c(rep("Oarfish (ONT)", 69 * 3), rep("Oarfish (PB)", 69 * 3), rep("Isoquant (ONT)", 69 * 3), rep("Isoquant (PB)", 69 * 3), rep("Bambu (ONT)", 69 * 3), rep("Bambu (PB)", 69 * 3), rep("Kallisto (ONT)", 69 * 3), rep("Kallisto (PB)", 69 * 3), rep("Miniquant (ONT)", 69 * 3), rep("Miniquant (PB)", 69 * 3), rep("Isosceles (ONT)", 69 * 3), rep("Isosceles (PB)", 69 * 3))
  annotation_frame_relative$fc <- unname(unlist(lapply(prepped_counts, function(x) c(x[, 1] - x[, 4], x[, 2] - x[, 5], x[, 3] - x[, 6]))))
  annotation_frame_relative$molarity_fold_change_factor <- log(annotation_frame_relative$molarity_fold_change, base = 2)
  annotation_frame_relative$tech <- factor(annotation_frame_relative$tech, levels = c("Bambu (ONT)", "Isoquant (ONT)", "Isosceles (ONT)", "Kallisto (ONT)", "Miniquant (ONT)", "Oarfish (ONT)", "Bambu (PB)", "Isoquant (PB)", "Isosceles (PB)", "Kallisto (PB)", "Miniquant (PB)", "Oarfish (PB)"))
  annotation_frame_relative$length_group <- ifelse(annotation_frame_relative$len < 750, "Exonic length < 750", "Exonic length >= 750")
  relative <- annotation_frame_relative %>% ggplot(aes(x = len, y = fc, color = as.factor(molarity_fold_change_factor), shape = as.factor(sample))) + geom_point(size = 4, alpha = 0.7) + labs(x = "SIRV length", y = "Observed logFC", fill = "", color = "Theoretical logFC", shape = "Replicate") + ggpubfigs::theme_big_simple() + ggplot2::scale_color_manual(values = ggpubfigs::friendly_pals$nickel_five) + geom_hline(yintercept = -6, col = "#648FFF", lty = 2, linewidth = 1) + geom_hline(yintercept = -1, 
      col = "#FE6100", lty = 2, linewidth = 1) + geom_hline(yintercept = 0, col = "#785EF0", lty = 2, linewidth = 1) + geom_hline(yintercept = 4, col = "#FFB000", lty = 2, linewidth = 1) + geom_vline(xintercept = 1250, lty = 2) + stat_cor(aes(x = as.numeric(molarity_fold_change_factor), color = NULL, y = fc, label = ..r.label..), method = "pearson", show.legend = FALSE, label.x.npc = c(0.025), label.y.npc = seq.int(1, 0.9875, length.out = 3), size = 4) + facet_wrap(~tech, nrow = 2) + theme(axis.text.x = element_text(angle = 90, 
      hjust = 1, vjust = 0.5))

  if (!exists("relative")) {
    stop("Panel object not found: relative", call. = FALSE)
  }

  panel_plot <- get("relative")
  ggplot2::ggsave(output_pdf, panel_plot, dpi = 600, width = 12, height = 8)
  ggplot2::ggsave(output_svg, panel_plot, dpi = 600, width = 12, height = 8)
  return(0)
}

if (exists("snakemake")) {
  log <- file(snakemake@log[[1]], open = "wt")
  sink(log, type = "output")
  sink(log, type = "message")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(edgeR)
  library(ggplot2)
  library(ggpubfigs)
  library(ggpubr)
  library(glue)
  library(rcartocolor)
  library(tibble)
  library(tidyr)
  library(tximeta)
  library(DRIMSeq)
  library(stageR)
  library(stringr)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplotify)
  library(readxl)
  library(rtracklayer)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_S03/panel_B.pdf"
  output_svg <- "results/figure_S03/panel_B.svg"
}

plot_figure_s03_b(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

