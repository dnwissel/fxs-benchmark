log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

dump_salmon_illumina <- function(input_path, output_path) {
  cs <- catchSalmon(paths = c(input_path))
  data.frame(
    transcript_id = names(cs$counts[, 1]),
    count = unname(cs$counts[, 1] / cs$annotation$Overdispersion)
  ) %>% write_tsv(output_path)
  return(0)
}

suppressPackageStartupMessages({
  library(edgeR)
  library(readr)
  library(dplyr)
  library(SummarizedExperiment)
})

status <- dump_salmon_illumina(
  input_path = snakemake@input[[1]],
  output_path = snakemake@output[[1]]
)

sink()
sink()
