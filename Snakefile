from snakemake.utils import min_version


configfile: "config/config.yaml"


min_version(config["snakemake_min_version"])


singularity: f"docker://continuumio/miniconda3:{config['miniconda_version']}"


include: "workflow/rules/setup.smk"
include: "workflow/rules/quantify_single_cell.smk"
include: "workflow/rules/quantify_bulk.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/preprocess_single_cell.smk"
include: "workflow/rules/prepare_bulk.smk"
include: "workflow/rules/downsample_single_cell.smk"
include: "workflow/rules/downsample_bulk.smk"
include: "workflow/rules/discover.smk"


rule all:
    input:
        expand(
            "results/quantify_bulk_downsampled/format/{method}/1_15000000.0/{tech}/gencode/{file_type}_counts_formatted.tsv",
            method=[
                "run_bambu_lr",
                "run_isoquant_lr",
                "run_isosceles_lr",
                "run_kallisto_long_lr",
                "run_miniquant_lr",
                "run_oarfish_lr",
            ],
            tech=["ont", "pb"],
            file_type=["transcript", "gene"],
        ),
        expand(
            "results/quantify_bulk_downsampled/format/{method}/1_15000000.0/illumina/gencode/{file_type}_counts_formatted.tsv",
            method=["run_salmon_illumina_corrected", "run_salmon_illumina_tpm"],
            file_type=["transcript", "gene"],
        ),
        expand(
            "results/format_quantify_subsampled/gencode/{tech}/1_60000000.0/{method}/pseudobulk/{file_type}_counts_formatted.tsv",
            tech=["ont", "pb"],
            method=["bambu", "isosceles", "kallisto", "oarfish"],
            file_type=["transcript", "gene"],
        ),
        expand(
            "results/qc/qc_prepare_sampled_read_quality_frame/{tech}/1000000.0/{platform}/{sample}.tsv",
            tech=["ont", "pb"],
            platform=["bulk", "single_cell"],
            sample=config["sample_names"],
        ),
        expand(
            "results/qc/qc_prepare_sampled_read_quality_frame_short_read_{platform}/1000000.0/{platform}/{sample}.tsv",
            platform=["bulk", "single_cell"],
            sample=config["sample_names"],
        ),
        expand(
            "results/qc/qc_calculate_three_prime_bias_binned/{tech}/1000000.0/{platform}/{sample}/{bin}/three_prime_bias.geneBodyCoverage.txt",
            tech=["ont", "pb"],
            platform=["bulk", "single_cell"],
            sample=config["sample_names"],
            bin=["lt1kb", "1-2kb", "2-3kb", "3-6kb", "gt6kb"],
        ),
        expand(
            "results/qc/qc_calculate_three_prime_bias_binned_short_reads/1000000.0/{platform}/{sample}/{bin}/three_prime_bias.geneBodyCoverage.txt",
            platform=["bulk", "single_cell"],
            sample=config["sample_names"],
            bin=["lt1kb", "1-2kb", "2-3kb", "3-6kb", "gt6kb"],
        ),
        expand(
            "results/qc/qc_calculate_internal_priming/{tech}/1000000.0/{platform}/{sample}_summary.txt",
            tech=["ont", "pb"],
            platform=["bulk", "single_cell"],
            sample=config["sample_names"],
        ),
        expand(
            "results/evaluate_performance_tusco/{dataset_type}/{tech}/1_{read_number}.0/{sample}.stats",
            dataset_type=["bulk", "sc"],
            tech=["ont", "pb"],
            read_number=[
                "2500000",
                "5000000",
                "10000000",
                "15000000",
                "20000000",
                "30000000",
                "60000000",
            ],
            sample=config["sample_names"],
        ),
        expand(
            "results/evaluate_performance_real/{dataset_type}/{tech}/1_{read_number}.0/{sample}.stats",
            dataset_type=["bulk", "sc"],
            tech=["ont", "pb"],
            read_number=[
                "2500000",
                "5000000",
                "10000000",
                "15000000",
                "20000000",
                "30000000",
                "60000000",
            ],
            sample=config["sample_names"],
        ),
        expand(
            "results/quantify_bulk_downsampled/format/{method}/1_75000.0/{tech}/sirv/transcript_counts_formatted.tsv",
            method=[
                "run_bambu_lr",
                "run_isoquant_lr",
                "run_isosceles_lr",
                "run_kallisto_long_lr",
                "run_miniquant_lr",
                "run_oarfish_lr",
            ],
            tech=["ont", "pb"],
        ),
        expand(
            "results/quantify_bulk_downsampled/format/{method}/1_20000.0/{tech}/ercc/transcript_counts_formatted.tsv",
            method=[
                "run_bambu_lr",
                "run_isoquant_lr",
                "run_isosceles_lr",
                "run_kallisto_long_lr",
                "run_miniquant_lr",
                "run_oarfish_lr",
            ],
            tech=["ont", "pb"],
        ),
        expand(
            "results/benchmarks/quantify_bulk_downsampled/{method}/1_{read_number}.0/gencode/{sample}/{tech}/{sample}.log",
            method=[
                "run_bambu_lr",
                "run_isoquant_lr",
                "run_isosceles_lr",
                "run_kallisto_long_lr",
                "run_miniquant_lr",
                "run_oarfish_lr",
            ],
            read_number=["1000000", "2500000", "5000000", "10000000", "15000000"],
            tech=["ont", "pb"],
            sample=config["sample_names"],
        ),
        expand(
            "results/quantify_bulk_downsampled/run_oarfish_lr_bootstrapped/1_15000000.0/gencode/{sample}/{tech}/{sample}.quant",
            tech=["pb"],
            sample=config["sample_names"],
        ),
