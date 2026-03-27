filter_gencode <- function(input_path_transcriptome,
                         input_path_illumina,
                         input_path_pb,
                         input_path_ont,
                         output_path_input,
                         output_path_output,
                         prop_removed = 0.1,
                         seed = 42) {
  set.seed(seed)
  transcriptome <- import(input_path_transcriptome)
  illumina <- vroom::vroom(input_path_illumina)
  pb <- vroom::vroom(input_path_pb)
  ont <- vroom::vroom(input_path_ont)
  
  expressed_transcripts <- Reduce(
    intersect,
    list(
      illumina$transcript_id[which(apply(edgeR::cpm(illumina[, 3:5]), 1, function(x) all(x >= 1)))],
      pb$transcript_id[which(apply(edgeR::cpm(pb[, 3:5]), 1, function(x) all(x >= 1)))],
      ont$transcript_id[which(apply(edgeR::cpm(ont[, 3:5]), 1, function(x) all(x >= 1)))]
    )
  )
  
  n_to_sample <- floor(prop_removed * length(expressed_transcripts))
  
  removed_transcripts <- sample(
    expressed_transcripts,
    n_to_sample
  )
  
  gencode_filtered <- transcriptome[!(transcriptome$transcript_id %in% removed_transcripts)]
  novel_filtered <- transcriptome[(transcriptome$transcript_id %in% removed_transcripts)]
  
  export(gencode_filtered, output_path_input)
  export(novel_filtered, output_path_output)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")
suppressPackageStartupMessages({
  library(rtracklayer)
  library(vroom)
  library(edgeR)
})

status <- filter_gencode(
  input_path_transcriptome = snakemake@input[["transcriptome"]],
  input_path_illumina = snakemake@input[["illumina"]],
  input_path_pb = snakemake@input[["kinnex"]],
  input_path_ont = snakemake@input[["ont"]],
  output_path_input = snakemake@output[[1]],
  output_path_output = snakemake@output[[2]]
)

sink()
sink()
