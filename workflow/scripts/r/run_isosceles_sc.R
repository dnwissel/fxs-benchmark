isosceles_quantify <- function(input_path_gtf,
                               input_path_fa,
                               input_path_bam,
                               threads,
                               output_path) {
  transcript_data <- prepare_transcripts(
    gtf_file = input_path_gtf,
    genome_fasta_file = input_path_fa
  )

  names(input_path_bam) <- "sample"

  se_tcc <- bam_to_tcc(
    input_path_bam,
    transcript_data,
    run_mode = "strict",
    min_read_count = 1,
    min_relative_expression = 0,
    is_single_cell = TRUE,
    ncpu = threads,
    barcode_tag = "BC"
  )

  se_quant <- tcc_to_transcript(
    se_tcc,
    ncpu = threads,
    use_length_normalization = FALSE
  )
  
  write_rds(
    se_quant,
    output_path
  )
  
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(BiocGenerics)
  library(rtracklayer)
  library(Isosceles)
  library(readr)
})

status <- isosceles_quantify(
  input_path_gtf = snakemake@input[["gencode_transcriptome"]],
  input_path_fa = snakemake@input[["gencode_genome"]],
  input_path_bam = snakemake@input[["reads"]],
  threads = snakemake@threads[[1]],
  output_path = snakemake@output[[1]]
)



sink()
sink()
