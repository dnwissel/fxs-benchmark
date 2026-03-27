import sys

import pandas as pd
import pysam
from tqdm import tqdm

with open(snakemake.log[0], "w") as f:
    sys.stderr = f
    sys.stdout = f

    for module in sys.modules.values():
        if hasattr(module, "__version__"):
            print(module.__name__, module.__version__)
        else:
            print(module.__name__)

    samfile = pysam.AlignmentFile(snakemake.input[0], "rb")

    read1 = None
    read2 = None

    read_name = []
    tlen = []

    for read in tqdm(
        samfile.fetch(until_eof=True), desc="Processing transcriptome reads"
    ):
        if (
            not read.is_paired
            or read.mate_is_unmapped
            or read.is_duplicate
            or read.is_secondary
            or read.is_supplementary
        ):
            continue

        if read.is_read2:
            read2 = read
        else:
            read1 = read
            read2 = None
            continue

        if (
            read1 is not None
            and read2 is not None
            and read1.query_name == read2.query_name
        ):
            read_name.append(read1.query_name)
            tlen.append(read1.template_length)

    pd.DataFrame({"qname": read_name, "tlen": tlen}).to_csv(
        snakemake.output[0], sep="\t", index=False
    )
