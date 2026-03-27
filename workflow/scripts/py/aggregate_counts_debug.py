import argparse
import sys
from functools import reduce
from pathlib import Path
from typing import List

import pandas as pd
from typeguard import typechecked


def _parse_gtf_annotation(gtf_path: Path) -> pd.DataFrame:
    try:
        gtf_df = pd.read_csv(
            gtf_path,
            sep="\t",
            header=None,
            comment="#",
            names=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute",
            ],
        )
    except FileNotFoundError:
        print(f"Error: GTF annotation file not found at {gtf_path}", file=sys.stderr)
        sys.exit(1)

    gtf_df = gtf_df[gtf_df["feature"] == "transcript"].copy()

    gtf_df["transcript_id"] = gtf_df["attribute"].str.extract(
        r'transcript_id "([^"]+)"'
    )
    gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')
    gtf_df.dropna(subset=["transcript_id", "gene_id"], inplace=True)
    gene_transcript_map = (
        gtf_df[["gene_id", "transcript_id"]].drop_duplicates().reset_index(drop=True)
    )

    gene_transcript_map["transcript_id"] = gene_transcript_map[
        "transcript_id"
    ].str.strip()

    return gene_transcript_map


@typechecked
def aggregate_counts(
    sample_paths: List[str],
    sample_names: List[str],
    gtf_annotation_path: str,
    output_path_transcript: str,
    output_path_gene: str,
    transcript_id_col_ix: int,
    count_col_ix: int,
) -> int:
    gtf_path = Path(gtf_annotation_path)
    output_transcript_path = Path(output_path_transcript)
    output_gene_path = Path(output_path_gene)
    output_transcript_path.parent.mkdir(parents=True, exist_ok=True)
    output_gene_path.parent.mkdir(parents=True, exist_ok=True)
    gene_transcript_map = _parse_gtf_annotation(gtf_path)
    dfs = []
    transcript_col_zero_ix = transcript_id_col_ix - 1
    count_col_zero_ix = count_col_ix - 1

    for sample_path_str, sample_name in zip(sample_paths, sample_names):
        sample_path = Path(sample_path_str)
        try:
            df = pd.read_csv(sample_path, sep="\t")
            df_subset = df.iloc[:, [transcript_col_zero_ix, count_col_zero_ix]].copy()
            df_subset.columns = ["transcript_id", sample_name]
            df_subset["transcript_id"] = (
                df_subset["transcript_id"].str.split(" ").str[0]
            )

            dfs.append(df_subset)
        except (FileNotFoundError, IndexError) as e:
            return 1

    merged_counts_df = reduce(
        lambda left, right: pd.merge(left, right, on="transcript_id", how="outer"), dfs
    )

    final_df = pd.merge(
        gene_transcript_map,
        merged_counts_df,
        on="transcript_id",
        how="left",
    )

    final_df[sample_names] = final_df[sample_names].fillna(0)
    final_df = final_df[["gene_id", "transcript_id"] + sample_names]
    final_df = final_df.sort_values(["gene_id", "transcript_id"]).reset_index(drop=True)
    final_df.to_csv(output_transcript_path, sep="\t", header=True, index=False)
    gene_df = final_df.groupby("gene_id")[sample_names].sum().reset_index()
    gene_df.to_csv(output_gene_path, sep="\t", header=True, index=False)

    return 0


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--sample_paths",
        nargs="*",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--sample_names",
        nargs="*",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--gtf_annotation_path",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--output_path_transcript",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--output_path_gene",
        required=True,
        type=str,
    )
    parser.add_argument(
        "--transcript_id_col_ix",
        required=True,
        type=int,
    )
    parser.add_argument(
        "--count_col_ix",
        required=True,
        type=int,
    )

    args = parser.parse_args()

    if len(args.sample_paths) != len(args.sample_names):
        print(
            "Error: The number of --sample_paths must equal the number of --sample_names.",
            file=sys.stderr,
        )
        sys.exit(1)

    aggregate_counts(
        sample_paths=args.sample_paths,
        sample_names=args.sample_names,
        gtf_annotation_path=args.gtf_annotation_path,
        output_path_transcript=args.output_path_transcript,
        output_path_gene=args.output_path_gene,
        transcript_id_col_ix=args.transcript_id_col_ix,
        count_col_ix=args.count_col_ix,
    )


if __name__ == "__main__":
    main()
