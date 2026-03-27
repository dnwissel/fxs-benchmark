create_illumina_quality_frame <- function(
    first_lengths,
    second_lengths,
    first_quality,
    second_quality,
    alignment,
    output_path) {

  alignment <- data.table::fread(
    alignment
  )

  selected_names <- alignment$qname

  length_first <- data.table::fread(
    first_lengths
  )

  length_second <- data.table::fread(
    second_lengths
  )

  quality_first <- data.table::fread(
    first_quality
  )

  quality_second <- data.table::fread(
    second_quality
  )
  
  length_first <- length_first[match(selected_names, length_first$V1), ]
  length_second <- length_second[match(selected_names, length_second$V1), ]

  quality_first <- quality_first[match(selected_names, quality_first$V1), ]
  quality_second <- quality_second[match(selected_names, quality_second$V1), ]
  
  alignment <- alignment[match(selected_names, alignment$qname), ]

  data.frame(
    qname = selected_names,
    length = length_first$V2 + length_second$V2,
    quality = apply(cbind(quality_first$V2, quality_second$V2), 1, mean),
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

status <- create_illumina_quality_frame(
  first_lengths = snakemake@input[["first_lengths"]],
  second_lengths = snakemake@input[["second_lengths"]],
  first_quality = snakemake@input[["first_quality"]],
  second_quality = snakemake@input[["second_quality"]],
  alignment = snakemake@input[["alignment"]],
  output_path = snakemake@output[[1]]
)

sink()
sink()
