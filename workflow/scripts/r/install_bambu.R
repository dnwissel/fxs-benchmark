install_bambu <- function(version, outfile) {
  renv::install("bioc::GenomeInfoDb")
  version <- "c5e923dd13249ee5c7563f2e378c795024c06ee7"
  renv::install(paste0("GoekeLab/bambu@", version), type="binary")
  # https://stackoverflow.com/questions/23922497/create-a-touch-file-on-unix
  write.table(data.frame(), file=outfile, col.names=FALSE)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages(
    {
        library(renv)
        library(BiocManager)
    }
)


status <- install_bambu(
    version = snakemake@params[["version"]],
    outfile = snakemake@output[[1]]
)

sink()
sink()
