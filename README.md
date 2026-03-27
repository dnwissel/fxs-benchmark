# Benchmarking long-read RNA-seq across modalities, methods, and sequencing depth in iNeurons

## Abstract

Long-read RNA sequencing (lrRNA-seq) provides advantages for transcript discovery and quantification through the sequencing of full-length transcripts. Although recent benchmarks have evaluated long-read technologies and quantification tools, to the best of our knowledge, no study to date has jointly compared sequencing technology, quantification choice, and depth across both bulk and single-cell platforms. Here, we generate a matched dataset using NGN2-induced neurons derived from Fragile X syndrome and isogenic rescue lines, profiled with bulk and single-cell Illumina, Oxford Nanopore Technologies (ONT), and Pacific Biosciences (PB) Kinnex technologies. All platforms and technologies capture the expected FMR1 reactivation signal. We find that PB bulk under-detects and under-quantifies short transcripts (less than 1.25 kb), ONT bulk under-detects and under-quantifies long transcripts (greater than 5 kb), and single-cell long-read technologies a large number of single-cell specific transcripts associated with truncations. Across six bulk and four single-cell long-read quantification tools, Isosceles, Miniquant, and Oarfish provide the best compromise between computational efficiency and performance in terms of quantification accuracy as measured by spike-ins, comparisons to Illumina, and on consensus based downstream tasks such as differential transcript expression (DTE). Depth-equivalency analyses reveal that PB single-cell sequencing requires approximately three- to four-fold greater depth than bulk to reach comparable power for transcript discovery and differential transcript expression. Our work aims to offer practical guidance for study design, including the choice of technology, sequencing depth, and quantification method. In addition, we hope our data may serve a reference dataset to evaluate emerging long-read transcriptomic protocols and methods as well as more closely investigate FMR1 biology.

## Reproducibility

To reproduce our results, you need to have `snakemake>=9.11.2`, `apptainer` and a recent `conda` (or alternative frontend, such as `mamba`) version installed. Starting from the raw data you can easily reproduce all of our analyses by running Snakemake:

```sh
snakemake --use-conda --use-apptainer --cores 12
```

## Data availability

### Raw data

Raw data is open access. Raw data of all single-cell and bulk samples are openly available from ArrayExpress (single-cell accession: E-MTAB-16805; bulk accession: E-MTAB-16791).

### Intermediate data (count tables for all downsampled methods, and for our recommended full-depth methods)

Intermediate data are available from [Zenodo](https://doi.org/10.5281/zenodo.19209644).

## Contact

In case of any questions, please reach out to [David](mailto:dwissel@ethz.ch) or open an issue in this repo.
