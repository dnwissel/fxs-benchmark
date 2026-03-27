import re
import sys

import pandas as pd
import pysam
from tqdm import tqdm

with open(snakemake.log[0], "w") as f:
    sys.stderr = f
    sys.stdout = f

    # Log loaded modules
    for module in sys.modules.values():
        if hasattr(module, "__version__"):
            print(module.__name__, module.__version__)
        else:
            print(module.__name__)

    samfile = pysam.AlignmentFile(snakemake.input[0], "rb")

    junction_pattern = re.compile(snakemake.params["junction_regex"])

    read1 = None
    read2 = None

    read_name = []
    edit_pct = []
    splice_junctions = []

    for read in tqdm(samfile.fetch(until_eof=True), desc="Processing paired reads"):
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
            try:
                nm1 = read1.get_tag("NM")
                nm2 = read2.get_tag("NM")
            except KeyError:
                continue

            total_len = read1.query_alignment_length + read2.query_alignment_length
            edit = (nm1 + nm2) / total_len

            # Union of junctions in both reads
            junctions = set(
                junction_pattern.findall(read1.cigarstring)
                + junction_pattern.findall(read2.cigarstring)
            )

            read_name.append(read1.query_name)
            edit_pct.append(edit)
            splice_junctions.append(len(junctions))

    pd.DataFrame(
        {"qname": read_name, "edit_distance": edit_pct, "n_junctions": splice_junctions}
    ).to_csv(snakemake.output[0], sep="\t", index=False)
