filter_tusco <- function(input_path_transcriptome,
                      input_path_tusco,
                      output_path_transcriptome,
                      output_path_tusco) {
  transcriptome <- import(input_path_transcriptome)
  tusco <- vroom::vroom(input_path_tusco, col_names = FALSE)
  tusco_transcriptome <- transcriptome[
    sapply(strsplit(transcriptome$gene_id, "\\."), function(x) x[[1]]) %in%
      tusco$X1
  ]
  transcriptome_without_tusco <- transcriptome[
    !(sapply(strsplit(transcriptome$gene_id, "\\."), function(x) x[[1]]) %in%
      tusco$X1)
  ]

  export(tusco_transcriptome, output_path_tusco)
  export(transcriptome_without_tusco, output_path_transcriptome)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(vroom)
})

status <- filter_tusco(
  input_path_transcriptome = snakemake@input[["transcriptome"]],
  input_path_tusco = snakemake@input[["tusco"]],
  output_path_transcriptome = snakemake@output[["transcriptome_without_tusco"]],
  output_path_tusco = snakemake@output[["tusco"]]
)

sink()
sink()
