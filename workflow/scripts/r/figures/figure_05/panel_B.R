plot_figure_05_b <- function(output_pdf, output_svg) {

  parse_gffcompare_f1 <- function(filepath) {
      if (!file.exists(filepath)) 
          return(NA)
      lines <- readLines(filepath)
      t_line <- grep("^Intron chain level:", lines, value = TRUE)
      if (length(t_line) == 0) 
          return(NA)
      stats <- as.numeric(unlist(regmatches(t_line, gregexpr("[0-9.]+", t_line))))
      if ((stats[1] + stats[2]) == 0) 
          return(0)
      return(2 * ((stats[1] * stats[2])/(stats[1] + stats[2])))
  }
  replicates <- c("E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3")
  bulk_depths <- c("1_2500000.0", "1_5000000.0", "1_10000000.0", "1_15000000.0")
  pseudo_depths <- c("1_10000000.0", "1_20000000.0", "1_30000000.0", "1_60000000.0")
  df_bulk <- crossing(Technology = c("ONT", "PB"), Dataset_Type = "bulk", Depth_Label = bulk_depths, Replicate = replicates)
  df_pseudo <- crossing(Technology = c("ONT", "PB"), Dataset_Type = "sc", Depth_Label = pseudo_depths, Replicate = replicates)
  df <- bind_rows(df_bulk, df_pseudo) %>% mutate(Numeric_Depth = as.numeric(gsub("^1_|\\.0$", "", Depth_Label)), path = paste0("/Volumes/Backup/20251120/evaluate_performance_tusco/", tolower(Dataset_Type), "/", tolower(Technology), "/", Depth_Label, "/", Replicate, ".stats"))
  df <- df %>% mutate(Base_Score = case_when(Dataset_Type == "bulk" & Technology == "ONT" ~ 70 + (Numeric_Depth/1e+06), Dataset_Type == "bulk" & Technology == "PB" ~ 75 + (Numeric_Depth/1e+06), Dataset_Type == "sc" & Technology == "ONT" ~ 55 + (Numeric_Depth/4e+07) * 25, TRUE ~ 60 + (Numeric_Depth/4e+07) * 25), F1_Score = Base_Score + rnorm(n(), mean = 0, sd = 1))
  df$F1_Score <- sapply(df$path, parse_gffcompare_f1)
  df_summary <- df %>% filter(!is.na(F1_Score)) %>% mutate(Dataset_Type = ifelse(Dataset_Type == "bulk", "Bulk", "Pseudobulk")) %>% group_by(Technology, Dataset_Type, Numeric_Depth) %>% summarise(Mean_F1 = mean(F1_Score), SD_F1 = sd(F1_Score), .groups = "drop")
  bulk_targets <- df_summary %>% filter(Dataset_Type == "Bulk") %>% dplyr::select(Technology, Target_Depth = Numeric_Depth, Target_F1 = Mean_F1)
  crossing_points <- bulk_targets %>% rowwise() %>% mutate(sc_equivalent_depth = get_sc_equivalent(Technology, Target_F1)) %>% ungroup() %>% filter(!is.na(sc_equivalent_depth))
  depth_discovery <- ggplot(df_summary, aes(x = Numeric_Depth, y = Mean_F1, color = Dataset_Type, fill = Dataset_Type)) + geom_ribbon(aes(ymin = Mean_F1 - SD_F1, ymax = Mean_F1 + SD_F1), alpha = 0.1, color = NA) + geom_line(linewidth = 1.5) + geom_point(size = 3) + geom_segment(data = crossing_points, aes(x = sc_equivalent_depth, xend = sc_equivalent_depth, y = 50, yend = Target_F1, group = Technology), color = "grey30", linetype = "dashed", alpha = 0.6, linewidth = 0.8, inherit.aes = FALSE) + geom_segment(data = crossing_points, 
      aes(x = Target_Depth, xend = sc_equivalent_depth, y = Target_F1, yend = Target_F1, group = Technology), color = "grey30", linetype = "dotted", alpha = 0.6, linewidth = 0.8, inherit.aes = FALSE) + geom_text(data = crossing_points, aes(x = sc_equivalent_depth, y = 49, label = paste0(round(sc_equivalent_depth/1e+06, 1), "M")), color = "black", inherit.aes = FALSE, angle = 0, hjust = 0.5, vjust = 0, size = 5, fontface = "bold") + facet_wrap(~Technology) + scale_x_continuous(labels = label_number(scale = 1e-06, 
      suffix = "M"), breaks = pretty_breaks(n = 5)) + scale_y_continuous(limits = c(49, 75)) + scale_color_manual(values = c(Bulk = "#E7B800", Pseudobulk = "#2E9FDF")) + scale_fill_manual(values = c(Bulk = "#E7B800", Pseudobulk = "#2E9FDF")) + theme_bw(base_size = 16) + labs(title = "Depth equivalency: transcript discovery (TUSCO)", y = "F1 Score (Intron Chain)", x = "Sequencing depth (million reads)", color = "Dataset", fill = "Dataset") + theme(legend.position = "bottom", strip.background = element_rect(fill = "gray95"), 
      strip.text = element_text(face = "bold"))

  if (!exists("depth_discovery")) {
    stop("Panel object not found: depth_discovery", call. = FALSE)
  }

  panel_plot <- get("depth_discovery")
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
  library(scales)
  library(edgeR)
  library(tidyr)
  library(ggplot2)
  library(dplyr)
  library(sessioninfo)
})

sessioninfo::session_info()

if (exists("snakemake")) {
  output_pdf <- snakemake@output[[1]]
  output_svg <- if (length(snakemake@output) >= 2) snakemake@output[[2]] else sub("\\.pdf$", ".svg", output_pdf)
} else {
  output_pdf <- "results/figure_05/panel_B.pdf"
  output_svg <- "results/figure_05/panel_B.svg"
}

plot_figure_05_b(output_pdf = output_pdf, output_svg = output_svg)

if (exists("snakemake")) {
  sink(type = "message")
  sink(type = "output")
  close(log)
}

