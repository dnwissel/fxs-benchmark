plot_figure_03_c <- function(output_pdf, output_svg) {

  config <- list(log_expr_pseudo = 0.1)
  cols <- carto_pal(7, "Safe")
  num_quant_methods_bulk <- 7
  cols <- carto_pal(num_quant_methods_bulk, "Safe")
  cols <- carto_pal(7, "Safe")
  plt_frame <- data.frame()
  for (depth in c("1000000.0", "2500000.0", "5000000.0", "10000000.0", "15000000.0")) {
      for (method in c("run_oarfish_lr", "run_miniquant_lr", "run_isoquant_lr", "run_isosceles_lr", "run_kallisto_long_lr", "run_bambu_lr")) {
          for (type in c("ont", "pb")) {
              for (rep in c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")) {
                  plt_frame <- rbind(plt_frame, data.frame(method = method, type = type, rep = rep, depth = as.numeric(depth), time = vroom::vroom(paste0("results/benchmarks/quantify_bulk_downsampled/", method, "/", "1_", depth, "/gencode/", rep, "/", type, "/", rep, ".log"))$s, mem = vroom::vroom(paste0("results/benchmarks/quantify_bulk_downsampled/", method, "/", "1_", depth, "/gencode/", rep, "/", type, "/", rep, ".log"))$max_uss))
              }
          }
      }
  }
  plt_frame$type <- ifelse(plt_frame$type == "ont", "ONT", "PB")
  plt_frame$method <- ifelse(plt_frame$method == "run_bambu_lr", "Bambu", ifelse(plt_frame$method == "run_isoquant_lr", "Isoquant", ifelse(plt_frame$method == "run_kallisto_long_lr", "Kallisto", ifelse(plt_frame$method == "run_oarfish_lr", "Oarfish", ifelse(plt_frame$method == "run_miniquant_lr", "Miniquant", ifelse(plt_frame$method == "run_isosceles_lr", "Isosceles", "Salmon"))))))
  plt_frame$method <- factor(plt_frame$method, levels = c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish", "Salmon"))
  plt_frame$type <- factor(plt_frame$type, levels = c("ONT", "PB"))
  plt_frame$depth <- factor(sprintf("%.1e", plt_frame$depth), levels = c("1.0e+06", "2.5e+06", "5.0e+06", "1.0e+07", "1.5e+07"))
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
  au_prcs <- unlist(lapply(list(ont_oarfish_counts, kinnex_oarfish_counts, ont_isoquant_counts, kinnex_isoquant_counts, ont_bambu_counts, kinnex_bambu_counts, ont_kallisto_counts, kinnex_kallisto_counts, ont_miniquant_counts, kinnex_miniquant_counts, ont_isosceles_counts, kinnex_isosceles_counts), function(counts) {
      lapply(c(0.01, 0.05, 0.1), function(cutoff) {
          dge <- edgeR::DGEList(counts = counts[, -(1:2)], genes = counts[, 2], group = c(0, 0, 0, 1, 1, 1))
          grp <- sapply(colnames(dge), function(x) strsplit(x, "-")[[1]][1])
          mm <- model.matrix(~grp)
          dge <- edgeR::estimateDisp(dge, mm)
          fit <- edgeR::glmQLFit(dge, mm)
          mc <- limma::makeContrasts(grpisoB11, levels = colnames(fit$coefficients))
          qlf <- edgeR::glmQLFTest(fit, contrast = mc)
          tt <- data.frame(edgeR::topTags(qlf, n = Inf)) %>% arrange(desc(transcript_id))
          truth <- annotation_frame_relative %>% dplyr::select(transcript_id, molarity_fold_change_factor) %>% distinct() %>% arrange(desc(transcript_id)) %>% mutate() %>% pull(molarity_fold_change_factor) %>% sign()
          prediction <- as.vector(decideTests(qlf, p.value = cutoff))
          positive_mask <- abs(truth) > 0
          tpr <- sum(truth[positive_mask] == (-1 * prediction[positive_mask]))/sum(positive_mask)
          fdr <- (sum(abs(prediction[!positive_mask])) + sum(abs(prediction[positive_mask][truth[positive_mask] == (prediction[positive_mask])])))/sum(abs(prediction))
          list(tpr, fdr)
      })
  }))
  plt_frame <- data.frame(tpr = au_prcs[seq.int(1, length(au_prcs), 2)], fdr = au_prcs[seq.int(2, length(au_prcs), 2)], method = rep(c("Oarfish (ONT)", "Oarfish (PB)", "Isoquant (ONT)", "Isoquant (PB)", "Bambu (ONT)", "Bambu (PB)", "Kallisto (ONT)", "Kallisto (PB)", "Miniquant (ONT)", "Miniquant (PB)", "Isosceles (ONT)", "Isosceles (PB)"), each = 3), rep = rep(rep(1:5, each = 3), each = 12), cutoff = rep(rep(c(0.01, 0.05, 0.1), 12), 5))
  plt_frame$method <- factor(plt_frame$method, levels = c("Bambu (ONT)", "Bambu (PB)", "Isoquant (ONT)", "Isoquant (PB)", "Isosceles (ONT)", "Isosceles (PB)", "Kallisto (ONT)", "Kallisto (PB)", "Miniquant (ONT)", "Miniquant (PB)", "Oarfish (ONT)", "Oarfish (PB)"))
  plt_frame$tech <- factor(ifelse(grepl("PB", plt_frame$method), "PB", "ONT"), levels = c("ONT", "PB"))
  plt_frame$method <- factor(ifelse(grepl("Bambu", plt_frame$method), "Bambu", ifelse(grepl("Kallisto", plt_frame$method), "Kallisto", ifelse(grepl("Miniquant", plt_frame$method), "Miniquant", ifelse(grepl("Oarfish", plt_frame$method), "Oarfish", ifelse(grepl("Isosceles", plt_frame$method), "Isosceles", "Isoquant"))))), levels = c("Bambu", "Isoquant", "Isosceles", "Kallisto", "Miniquant", "Oarfish"))
  dte_sirvs <- plt_frame %>% group_by(cutoff, method, tech) %>% summarise(tpr = mean(tpr), fdr = mean(fdr)) %>% ggplot(aes(x = fdr, y = tpr, color = method, group = method)) + geom_point(size = 4) + geom_path(linewidth = 1, lty = 2) + geom_vline(xintercept = 0.01, lty = 2) + geom_vline(xintercept = 0.05, lty = 2) + geom_vline(xintercept = 0.1, lty = 2) + theme_big_simple() + labs(fill = "", y = "TPR (SIRVs)", x = "FDR (SIRVs)", color = "") + scale_color_manual(values = rep(cols[1:6], each = 1)) + facet_wrap(~tech, 
      nrow = 2) + guides(colour = guide_legend(nrow = 2))

  if (!exists("dte_sirvs")) {
    stop("Panel object not found: dte_sirvs", call. = FALSE)
  }

  panel_plot <- get("dte_sirvs")
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
  output_pdf <- "results/figure_03/panel_C.pdf"
  output_svg <- "results/figure_03/panel_C.svg"
}

plot_figure_03_c(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

