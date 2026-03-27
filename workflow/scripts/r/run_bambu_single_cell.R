run_bambu <- function(reads, annotations, genome, discovery, output_path, ncore, NDR, quantification, lowMemory = TRUE) {
  bambuAnnotations <- prepareAnnotations(annotations)
  readClassFile <- bambu(
    reads = reads,
    annotations = bambuAnnotations,
    genome = genome,
    ncore = ncore,
    discovery = FALSE,
    quant = FALSE,
    demultiplexed = TRUE,
    verbose = TRUE,
    assignDist = FALSE,
    lowMemory = lowMemory,
    yieldSize = 10000000,
    cleanReads = TRUE,
    dedupUMI = TRUE
  )
  quantData <- bambu(
    reads = readClassFile, annotations = bambuAnnotations, genome = genome,
    ncore = ncore, discovery = FALSE, quant = FALSE, demultiplexed = TRUE, verbose = FALSE, opt.em = list(degradationBias = FALSE), assignDist = TRUE, spatial = NULL
  )
  bambu_result <- bambu(
    reads = NULL,
    annotations = bambuAnnotations,
    genome = genome,
    quantData = quantData,
    assignDist = FALSE,
    ncore = ncore,
    discovery = FALSE,
    quant = TRUE,
    demultiplexed = TRUE,
    verbose = FALSE,
    opt.em = list(degradationBias = FALSE),
    clusters = NULL
  )
  saveRDS(
	          bambu_result,
		      output_path
  )

}

suppressPackageStartupMessages(library(bambu))
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")


run_bambu(
  reads = snakemake@input[["reads"]],
  annotations = snakemake@input[["gencode_transcriptome"]],
  genome = snakemake@input[["gencode_genome"]],
  discovery = FALSE,
  output_path = snakemake@output[[1]],
  ncore = snakemake@threads,
  lowMemory = TRUE,
  NDR = NULL,
  quantification = TRUE
)
