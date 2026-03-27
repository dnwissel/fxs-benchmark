create_kinnex_quality_frame <- function(
    lengths,
    quality,
    alignment,
    output_path) {

  alignment <- data.table::fread(
    alignment
  )

  selected_names <- alignment$qname

  lengths <- data.table::fread(
    lengths
  )
  quality <- data.table::fread(
    quality
  )

  alignment <- alignment[match(selected_names, alignment$qname), ]
  lengths <- lengths[match(selected_names, lengths$V1), ]
  quality <- quality[match(selected_names, quality$V1), ]

  data.frame(
    qname = selected_names,
    length = lengths$V2,
    quality = quality$V2,
    n_junctions = alignment$n_junctions,
    edit_distance = alignment$edit_distance,
    snv_rate = alignment$snv_rate,
    indel_rate = alignment$indel_rate
  ) %>% write_tsv(output_path)
}



log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
})

status <- create_kinnex_quality_frame(
  lengths = snakemake@input[["lengths"]],
  quality = snakemake@input[["quality"]],
  alignment = snakemake@input[["alignment"]],
  output_path = snakemake@output[[1]]
)

sink()
sink()
