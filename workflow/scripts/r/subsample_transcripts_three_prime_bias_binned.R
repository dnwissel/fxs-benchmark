log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

subsample_transcripts_three_prime_bias_binned <- function(
    transcriptome,
    less_one_kb,
    one_two_kb,
    two_three_kb,
    three_six_kb,
    greater_six_kb,
    n_subsample,
    seed) {
  set.seed(as.integer(seed))
  transcriptome <- import(transcriptome)
  print(head(
    data.frame(transcriptome[transcriptome$type == "exon"]) %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::summarise(length = sum(width))
  ))
  transcriptome_grouped <- data.frame(transcriptome[transcriptome$type == "exon"]) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(length = sum(width)) %>%
    dplyr::mutate(group = ifelse(
      length <= 1000,
      "<1kb",
      ifelse(
        length <= 2000,
        "1-2kb",
        ifelse(
          length <= 3000,
          "2-3kb",
          ifelse(
            length <= 6000,
            "3-6kb",
            ">6kb"
          )
        )
      )
    ))

  transcript_ids_one <- transcriptome_grouped %>%
    filter(group == "<1kb") %>%
    pull(transcript_id) %>%
    unique()
  transcript_one_two <- transcriptome_grouped %>%
    filter(group == "1-2kb") %>%
    pull(transcript_id) %>%
    unique()
  transcript_ids_two_three <- transcriptome_grouped %>%
    filter(group == "2-3kb") %>%
    pull(transcript_id) %>%
    unique()
  transcript_ids_three_six <- transcriptome_grouped %>%
    filter(group == "3-6kb") %>%
    pull(transcript_id) %>%
    unique()
  transcript_ids_six <- transcriptome_grouped %>%
    filter(group == ">6kb") %>%
    pull(transcript_id) %>%
    unique()

  sampled_transcripts_one <- sample(x = transcript_ids_one, replace = FALSE, size = n_subsample)
  sampled_transcripts_one_two <- sample(x = transcript_one_two, replace = FALSE, size = n_subsample)
  sampled_transcripts_two_three <- sample(x = transcript_ids_two_three, replace = FALSE, size = n_subsample)
  sampled_transcripts_three_six <- sample(x = transcript_ids_three_six, replace = FALSE, size = n_subsample)
  sampled_transcripts_six <- sample(x = transcript_ids_six, replace = FALSE, size = n_subsample)


  export(
    transcriptome[transcriptome$transcript_id %in% sampled_transcripts_one],
    less_one_kb
  )
  export(
    transcriptome[transcriptome$transcript_id %in% sampled_transcripts_one_two],
    one_two_kb
  )
  export(
    transcriptome[transcriptome$transcript_id %in% sampled_transcripts_two_three],
    two_three_kb
  )
  export(
    transcriptome[transcriptome$transcript_id %in% sampled_transcripts_three_six],
    three_six_kb
  )
  export(
    transcriptome[transcriptome$transcript_id %in% sampled_transcripts_six],
    greater_six_kb
  )

  return(0)
}

suppressPackageStartupMessages({
  library(rtracklayer)
  library(dplyr)
})

status <- subsample_transcripts_three_prime_bias_binned(
  transcriptome = snakemake@input[["transcriptome"]],
  less_one_kb = snakemake@output[["less_one_kb"]],
  one_two_kb = snakemake@output[["one_two_kb"]],
  two_three_kb = snakemake@output[["two_three_kb"]],
  three_six_kb = snakemake@output[["three_six_kb"]],
  greater_six_kb = snakemake@output[["greater_six_kb"]],
  n_subsample = snakemake@params[["n_subsample"]],
  seed = snakemake@params[["seed"]]
)

sessionInfo()

sink()
sink()
