run_bambu <- function(reads, annotations, genome, discovery, output_directory, ncore, NDR, quantification, lowMemory = TRUE) {
    bambuAnnotations <- prepareAnnotations(annotations)
    bambu_result <- bambu::bambu(reads = reads, annotations = bambuAnnotations,
                                 genome = genome, discovery = discovery,
                                 ncore = ncore, lowMemory = lowMemory, NDR = NDR, quant = quantification
    )
    bambu::writeBambuOutput(bambu_result, path = output_directory)
}

suppressPackageStartupMessages(library(bambu))
log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")


if (snakemake@wildcards[["data_type"]] == "gencode")  {
    run_bambu(
        reads = snakemake@input[["reads"]],
        annotations = snakemake@input[["gencode_transcriptome"]],
        genome = snakemake@input[["gencode_genome"]],
        discovery = snakemake@params[["discovery"]],
        output_directory = snakemake@params[["outdir"]],
        ncore = snakemake@threads,
        lowMemory = TRUE,
        NDR = snakemake@params[["NDR"]],
        quantification = snakemake@params[["quantification"]]
        )
} else {
    run_bambu(
        reads = snakemake@input[["reads"]],
        annotations = snakemake@input[["sirv_transcriptome"]],
        genome = snakemake@input[["sirv_genome"]],
        discovery = snakemake@params[["discovery"]],
        output_directory = snakemake@params[["outdir"]],
        ncore = snakemake@threads,
        lowMemory = TRUE,
        NDR = snakemake@params[["NDR"]],
        quantification = snakemake@params[["quantification"]]
        )
}
