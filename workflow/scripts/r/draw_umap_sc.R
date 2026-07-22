log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(bluster)
  library(SingleCellExperiment)
  library(scran)
  library(ggplot2)
  library(dplyr)
  library(scater)
})

ill_sce <- readRDS(snakemake@input[["sce"]])
plotdir <- snakemake@output[["plotdir"]]
meta_path <- snakemake@input[["md"]]
meta_full <- read.csv(meta_path, header = TRUE) 
meta <- meta_full %>% filter(method == "Illumina" & type == "sc")

# Leiden clustering and fates determination
set.seed(123)
ill_sce$leiden <- clusterCells(ill_sce,
                               use.dimred = "HARMONY",
                               BLUSPARAM = SNNGraphParam(
                                k = snakemake@params[["k"]],
                                cluster.fun = "leiden"
                                )
)

umap_ids <- plotUMAP(ill_sce, colour_by = "sample_name", point_size = 1.5) +
  scale_color_manual(values = c(
    "E3-1" = "#4352F0",
    "E3-2" = "#F04443",
    "E3-3" = "#F0CB43",
    "isoB11-1" = "#43F075",
    "isoB11-2" = "#58705F",
    "isoB11-3" = "#5F649B" 
  )) +
  labs(color = "Sample id") +
  theme_minimal() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        legend.position = c(0.95, 0.95),  
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 1), color = NA),
        legend.box.background = element_rect(color = "black", linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)
  )

# save umap ids
save_path <- file.path(plotdir, "umap_ids.pdf")
ggsave(filename = save_path,
       plot = umap_ids,
       width = 16,
       height = 12,
       units = "cm",
       dpi = 300,
       device = "pdf")
save_path <- file.path(plotdir, "umap_ids.svg")
ggsave(filename = save_path,
       plot = umap_ids,
       width = 16,
       height = 12,
       units = "cm",
       dpi = 300,
       device = "svg")

umap_clusters <- plotUMAP(ill_sce, colour_by = "leiden", point_size = 1.5) +
  scale_color_manual(values = c(
    "1" = "#4352F0",
    "2" = "#F04443",
    "3" = "#F0CB43",
    "4" = "#43F075",
    "5" = "#58705F",
    "6" = "#5F649B" 
  )) +
  labs(color = "Leiden clusters") +
  theme_minimal() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        legend.position = c(0.95, 0.95),  
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 1), color = NA),
        legend.box.background = element_rect(color = "black", linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)
  )

# save umap clusters
save_path <- file.path(plotdir, "umap_clusters.pdf")
ggsave(filename = save_path,
       plot = umap_clusters,
       width = 16,
       height = 12,
       units = "cm",
       dpi = 300,
       device = "pdf")
save_path <- file.path(plotdir, "umap_clusters.svg")
ggsave(filename = save_path,
       plot = umap_clusters,
       width = 16,
       height = 12,
       units = "cm",
       dpi = 300,
       device = "svg")

ill_sce$fates <- case_when(
  ill_sce$leiden == "1" ~ "Fate 3",
  ill_sce$leiden == "2" ~ "Fate 2",
  TRUE ~ "Other fate"
)
umap_fates <- plotUMAP(ill_sce, colour_by = "fates", point_size = 1.5) +
  scale_color_manual(values = c(
    "Fate 3" = "#4352F0",
    "Fate 2" = "#F04443",
    "Other fate" = "grey"
  )) +
  labs(color = "Neuronal fates") +
  theme_minimal() +
  theme(legend.text = element_text(size=12),
        legend.title = element_text(size=16),
        legend.position = c(0.95, 0.95),  
        legend.justification = c("right", "top"),
        legend.background = element_rect(fill = alpha("white", 1), color = NA),
        legend.box.background = element_rect(color = "black", linewidth = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 12)
  )

# save umap fates
save_path <- file.path(plotdir, "umap_fates.pdf")
ggsave(filename = save_path,
       plot = umap_fates,
       width = 16,
       height = 12,
       units = "cm",
       dpi = 300,
       device = "pdf")
save_path <- file.path(plotdir, "umap_fates.svg")
ggsave(filename = save_path,
       plot = umap_fates,
       width = 16,
       height = 12,
       units = "cm",
       dpi = 300,
       device = "svg")

### TODO: decide how to plot these
# reproduce Treutlein's paper markers
# NANOG
umap_nanog <- plotUMAP(ill_sce, colour_by = "ENSG00000111704", point_size = 1.5, point_alpha=0.5) +
  ggtitle("NANOG expression")
grid.arrange(umap_nanog, umap_fates, ncol=2)

# MAP2
umap_map2 <- plotUMAP(ill_sce, colour_by = "ENSG00000078018", point_size = 1.5, point_alpha=0.5) +
  ggtitle("MAP2 expression")
grid.arrange(umap_map2, umap_clusters, ncol=2)

# PRPH
umap_prph <- plotUMAP(ill_sce, colour_by = "ENSG00000135406", point_size = 1.5, point_alpha=0.5) +
  ggtitle("PRPH expression")
grid.arrange(umap_prph, umap_clusters, ncol=2)

# PHOX2B
# no PHOX2B in filtered genes

# POU4F1
umap_pou4f1 <- plotUMAP(ill_sce, colour_by = "ENSG00000152192", point_size = 1.5, point_alpha=0.5) +
  ggtitle("POU4F1 expression")
grid.arrange(umap_pou4f1, umap_clusters, ncol=2)

# LHX9
umap_lhx9 <- plotUMAP(ill_sce, colour_by = "ENSG00000143355", point_size = 1.5, point_alpha=0.5) +
  ggtitle("LHX9 expression")
grid.arrange(umap_lhx9, umap_clusters, ncol=2)

# GPM6A
umap_gpm6a <- plotUMAP(ill_sce, colour_by = "ENSG00000150625", point_size = 1.5, point_alpha=0.5) +
  ggtitle("GPM6A expression")
grid.arrange(umap_gpm6a, umap_clusters, ncol=2)

# Probably worth adding some general markers for excitatory neurons and GRIA2 (as Will suggested)

sink()
sink()