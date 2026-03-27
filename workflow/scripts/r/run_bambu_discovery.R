run_bambu <- function(reads, annotations, genome, discovery, output_directory, ncore, NDR, quantification, lowMemory = TRUE) {
    bambuAnnotations <- prepareAnnotations(annotations)
    bambu_result <- bambu::bambu(reads = reads, annotations = bambuAnnotations,
                                 genome = genome, discovery = discovery,
                                 ncore = ncore, lowMemory = lowMemory, NDR = NDR, quant = quantification
    )
    writeToGTF(rowRanges(bambu_result[[1]]), output_directory)
}

suppressPackageStartupMessages(library(bambu))
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

run_bambu(
    reads = snakemake@input[["reads"]],
    annotations = snakemake@input[["gencode_transcriptome"]],
    genome = snakemake@input[["gencode_genome"]],
    discovery = TRUE,
    output_directory = snakemake@output[[1]],
    ncore = snakemake@threads,
    lowMemory = TRUE,
    NDR = 0.05,
    quantification = FALSE
    )
