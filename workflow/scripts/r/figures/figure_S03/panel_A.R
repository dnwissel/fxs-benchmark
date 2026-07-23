plot_figure_s03_a <- function(output_pdf, output_svg) {

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
  annotation_frame_absolute <- data.frame(transcript_id = annot$`GenBank*`[2:93], mw_g_mol = as.numeric(annot[, 7][2:93]), stock_concentration = as.numeric(annot[, 5][2:93]), conc = as.numeric(annot[, 6][2:93]), conc_ng_ul = as.numeric(annot[, 8][2:93]), length = as.numeric(annot[, 9][2:93])) %>% arrange(desc(transcript_id))
  annotation_frame_absolute <- rbind(annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute)
  annotation_frame_absolute$sample <- rep(c("E3-1", "E3-2", "E3-3"), each = 92)
  annotation_frame_absolute <- rbind(annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute, annotation_frame_absolute)
  annotation_frame_absolute$tech <- c(rep("Bambu (ONT)", 92 * 3), rep("Bambu (PB)", 92 * 3), rep("Isoquant (ONT)", 92 * 3), rep("Isoquant (PB)", 92 * 3), rep("Isosceles (ONT)", 92 * 3), rep("Isosceles (PB)", 92 * 3), rep("Kallisto (ONT)", 92 * 3), rep("Kallisto (PB)", 92 * 3), rep("Miniquant (ONT)", 92 * 3), rep("Miniquant (PB)", 92 * 3), rep("Oarfish (ONT)", 92 * 3), rep("Oarfish (PB)", 92 * 3))
  annotation_frame_absolute$tech <- factor(annotation_frame_absolute$tech, levels = c("Bambu (ONT)", "Isoquant (ONT)", "Isosceles (ONT)", "Kallisto (ONT)", "Miniquant (ONT)", "Oarfish (ONT)", "Bambu (PB)", "Isoquant (PB)", "Isosceles (PB)", "Kallisto (PB)", "Miniquant (PB)", "Oarfish (PB)"))
  annotation_frame_absolute$log_cpm <- unlist(lapply(prepped_counts, function(x) {
      c(x[, 1], x[, 2], x[, 3])
  }))
  annotation_frame_absolute$rep <- sapply(strsplit(annotation_frame_absolute$sample, "\\-"), function(x) x[[2]])
  annotation_frame_absolute$day <- sapply(strsplit(annotation_frame_absolute$sample, "\\-"), function(x) x[[1]])
  annotation_frame_absolute_ont <- annotation_frame_absolute
  annotation_frame_absolute_ont$pred_concentration <- unlist(lapply(as.character(unique(annotation_frame_absolute_ont$tech)), function(local_tech) {
      lapply(c("1", "2", "3"), function(local_rep) {
          annotation_frame_absolute_local <- annotation_frame_absolute_ont %>% filter(rep == local_rep & tech == local_tech)
          reg <- lm(log(annotation_frame_absolute_local$stock_concentration + 1) ~ annotation_frame_absolute_local$log_cpm)
          predict(reg, annotation_frame_absolute_local)
      })
  }))
  annotation_frame_absolute_ont$length_cat <- factor(ifelse(annotation_frame_absolute_ont$length <= 1250, "<= 1250", "> 1250"), levels = c("<= 1250", "> 1250"))
  absolute <- annotation_frame_absolute_ont %>% ggplot(aes(y = log(conc + 1), x = pred_concentration, shape = as.factor(rep))) + geom_point(size = 3) + labs(x = "Predicted concentration (logCPM)", y = "Log concentration (amoles/\302\265l)", fill = "", shape = "Replicate") + theme_big_simple() + geom_smooth(aes(shape = NULL), method = "lm", se = FALSE, color = "lightblue") + stat_cor(aes(x = as.numeric(pred_concentration), y = log(stock_concentration + 1), color = NULL, label = ..r.label..), method = "pearson", 
      show.legend = FALSE, label.x.npc = c(0.025), label.y.npc = seq.int(1, 0.975 - 0.0125, length.out = 3), size = 4) + facet_grid(vars(length_cat), vars(tech)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

  if (!exists("absolute")) {
    stop("Panel object not found: absolute", call. = FALSE)
  }

  panel_plot <- get("absolute")
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
  output_pdf <- "results/figure_S03/panel_A.pdf"
  output_svg <- "results/figure_S03/panel_A.svg"
}

plot_figure_s03_a(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

