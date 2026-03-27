import re
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

    junction_pattern = re.compile(snakemake.params["junction_regex"])

    read_name = []
    edit_pct = []
    splice_junctions = []

    for read in tqdm(samfile.fetch(until_eof=True), desc="Processing reads"):
        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue

        try:
            nm = read.get_tag("NM")
            length = read.query_alignment_length
        except KeyError:
            continue

        read_name.append(read.query_name)
        edit_pct.append(nm / length)
        splice_junctions.append(len(junction_pattern.findall(read.cigarstring)))

    pd.DataFrame(
        {"qname": read_name, "edit_distance": edit_pct, "n_junctions": splice_junctions}
    ).to_csv(snakemake.output[0], sep="\t", index=False)
