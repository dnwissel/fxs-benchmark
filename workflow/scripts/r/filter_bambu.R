filter_bambu <- function(input_path_transcriptome,
                      input_path_tusco,
                      output_path) {
  transcriptome <- import(input_path_transcriptome)
  tusco <- import(input_path_tusco)
  overlaps <- findOverlaps(transcriptome, tusco)
  indices_to_keep <- unique(queryHits(overlaps))
  filtered_transcriptome <- transcriptome[indices_to_keep]
  export(filtered_transcriptome, output_path)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(vroom)
})

status <- filter_bambu(
  input_path_transcriptome = snakemake@input[["transcriptome"]],
  input_path_tusco = snakemake@input[["tusco"]],
  output_path = snakemake@output[[1]]
)

sink()
sink()
