filter_bambu_novel <- function(input_path_transcriptome,
                         output_path) {
  transcriptome <- import(input_path_transcriptome)
  filtered_transcriptome <- transcriptome[
    grepl("Bambu", transcriptome$transcript_id)
  ]
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

status <- filter_bambu_novel(
  input_path_transcriptome = snakemake@input[["transcriptome"]],
  output_path = snakemake@output[[1]]
)

sink()
sink()
