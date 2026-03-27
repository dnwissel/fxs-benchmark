# -*- coding: utf-8 -*-

import re
import sys

import pandas as pd
import pysam


def calculate_snv_indel_rates(read):
    if not read.has_tag("MD"):
        return None, None

    aligned_length = read.query_alignment_length
    if aligned_length == 0:
        return None, None

    md_tag = read.get_tag("MD")
    parts = re.split(r"(\d+)", md_tag)
    snv_count = 0
    for part in parts:
        if not part.isdigit() and part and not part.startswith("^"):
            snv_count += len(part)

    indel_length = 0
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == 1 or op == 2:
                indel_length += length

    snv_rate = snv_count / aligned_length
    indel_rate = indel_length / aligned_length

    return snv_rate, indel_rate


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    print("--- Package Versions ---")
    for module in sys.modules.values():
        if hasattr(module, "__version__"):
            print(f"{module.__name__}: {module.__version__}")
    print("------------------------\n")

    try:

        junction_regex = snakemake.params.get("junction_regex", "N")
        data_records = []

        with pysam.AlignmentFile(snakemake.input[0], "rb") as samfile:
            for i, read in enumerate(samfile):
                if (
                    read.is_duplicate
                    or read.is_secondary
                    or read.is_supplementary
                    or read.is_unmapped
                ):
                    continue

                try:
                    aligned_length = read.query_alignment_length
                    if aligned_length == 0:
                        continue

                    snv_rate, indel_rate = calculate_snv_indel_rates(read)
                    if snv_rate is None:
                        continue

                    edit_distance_pct = read.get_tag("NM") / aligned_length
                    n_junctions = len(re.findall(junction_regex, read.cigarstring))

                    data_records.append(
                        {
                            "qname": read.query_name,
                            "edit_distance": edit_distance_pct,
                            "n_junctions": n_junctions,
                            "snv_rate": snv_rate,
                            "indel_rate": indel_rate,
                        }
                    )

                except (KeyError, ZeroDivisionError):
                    continue

                if (i + 1) % 500000 == 0:
                    print(f"  ...processed {i+1} reads")

        df = pd.DataFrame.from_records(data_records)

        column_order = [
            "qname",
            "edit_distance",
            "n_junctions",
            "snv_rate",
            "indel_rate",
        ]
        df = df[column_order]

        df.to_csv(snakemake.output[0], sep="\t", index=False, float_format="%.6f")

    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
        sys.exit(1)
