# -*- coding: utf-8 -*-

import re
import sys
import pandas as pd
import pysam


def get_error_counts(read):
    if not read.has_tag('MD') or not read.has_tag('NM') or read.is_unmapped:
        return None

    aligned_length = read.query_alignment_length
    if aligned_length == 0:
        return None

    md_tag = read.get_tag('MD')
    parts = re.split(r'(\d+)', md_tag)
    snv_count = 0
    for part in parts:
        if not part.isdigit() and part and not part.startswith('^'):
            snv_count += len(part)

    indel_length = 0
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == 1 or op == 2:
                indel_length += length
    
    return snv_count, indel_length, aligned_length


with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    
    print("--- Package Versions ---")
    for module in sys.modules.values():
        if hasattr(module, "__version__"):
            print(f"{module.__name__}: {module.__version__}")
    print("------------------------\n")

    try:
        samfile = pysam.AlignmentFile(snakemake.input[0], "rb")

        records = []
        
        read1 = None
        read2 = None
        
        junction_regex = snakemake.params.get("junction_regex", "N")

        for i, current_read in enumerate(samfile):
            if (
                not current_read.is_paired
                or current_read.mate_is_unmapped
                or current_read.is_duplicate
                or current_read.is_secondary
                or current_read.is_supplementary
            ):
                continue

            if current_read.is_read2:
                read2 = current_read
            else:
                read1 = current_read
                read2 = None
                continue

            if (read1 is not None and read2 is not None and read1.query_name == read2.query_name):
                
                try:
                    r1_counts = get_error_counts(read1)
                    r2_counts = get_error_counts(read2)
                    
                    if r1_counts is None or r2_counts is None:
                        continue
                    
                    snv1, indel1, len1 = r1_counts
                    snv2, indel2, len2 = r2_counts
                    
                    total_aligned_length = len1 + len2
                    if total_aligned_length == 0:
                        continue
                    total_nm = read1.get_tag("NM") + read2.get_tag("NM")
                    combined_edit_distance = total_nm / total_aligned_length
                    
                    total_snvs = snv1 + snv2
                    combined_snv_rate = total_snvs / total_aligned_length
                    
                    total_indels = indel1 + indel2
                    combined_indel_rate = total_indels / total_aligned_length

                    junctions1 = re.findall(junction_regex, read1.cigarstring)
                    junctions2 = re.findall(junction_regex, read2.cigarstring)
                    unique_junctions_count = len(set(junctions1 + junctions2))

                    records.append({
                        "qname": read1.query_name,
                        "edit_distance": combined_edit_distance,
                        "n_junctions": unique_junctions_count,
                        "snv_rate": combined_snv_rate,
                        "indel_rate": combined_indel_rate,
                    })

                except (KeyError, ZeroDivisionError):
                    continue
                
                read1, read2 = None, None

            if (i + 1) % 1000000 == 0:
                 print(f"  ...processed {i+1} total reads")


        
        if not records:
             print("Warning: No valid read pairs were found to process.")
        else:
            df = pd.DataFrame.from_records(records)
            
            column_order = ['qname', 'edit_distance', 'n_junctions', 'snv_rate', 'indel_rate']
            df = df[column_order]

            df.to_csv(
                snakemake.output[0], 
                sep="\t", 
                index=False, 
                float_format="%.6f"
            )

    except Exception as e:
        sys.exit(1)
