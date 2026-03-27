install_isosceles <- function(version, outfile) {
  renv::install(paste0("Genentech/Isosceles@", version), type = "binary")
  # https://stackoverflow.com/questions/23922497/create-a-touch-file-on-unix
  write.table(data.frame(), file = outfile, col.names = FALSE)
  return(0)
}

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

suppressPackageStartupMessages({
  library(renv)
})


status <- install_isosceles(
  version = snakemake@params[["version"]],
  outfile = snakemake@output[[1]]
)

sink()
sink()
