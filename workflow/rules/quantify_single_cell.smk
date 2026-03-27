rule quantify_sc_run_oarfish:
    input:
        "results/preprocess/sort_{tech}_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    output:
        "results/quantify_sc_run_oarfish/{tech}/{subsample_number}_{number_to_sample}/{sample}.count.mtx",
    log:
        "logs/quantify_sc_run_oarfish/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    benchmark:
        repeat(
            "benchmarks/preprocess/align_{tech}_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.log",
            config["repetitions"],
        )
    conda:
        "../envs/standalone/oarfish.yaml"
    threads: config["quantify_threads"]
    params:
        bin_width=100,
        filter_group="no-filters",
        output_prefix="results/quantify_sc_run_oarfish/{tech}/{subsample_number}_{number_to_sample}/{sample}",
    shell:
        """
        conda list &> {log};
        oarfish --threads {threads} \
            --filter-group {params.filter_group} \
            --model-coverage \
            --bin-width {params.bin_width} \
            --alignments {input} \
            --single-cell \
            --output {params.output_prefix} &>> {log}
        """


rule quantify_sc_run_oarfish_bootstrapped:
    input:
        "results/preprocess/sort_{tech}_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    output:
        "results/quantify_sc_run_oarfish_bootstrapped/{tech}/{subsample_number}_{number_to_sample}/{sample}.quant",
    log:
        "logs/quantify_sc_run_oarfish_bootstrapped/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/oarfish.yaml"
    threads: config["quantify_threads"]
    params:
        bin_width=100,
        filter_group="no-filters",
        n_bootstraps=30,
        output_prefix="results/quantify_sc_run_oarfish_bootstrapped/{tech}/{subsample_number}_{number_to_sample}/{sample}",
    shell:
        """
        conda list &> {log};
        oarfish --threads {threads} \
            --filter-group {params.filter_group} \
            --model-coverage \
            --bin-width {params.bin_width} \
            --num-bootstraps {params.n_bootstraps} \
            --alignments {input} \
            --output {params.output_prefix} &>> {log}
        """


rule quantify_sc_run_bambu:
    input:
        "results/setup/install_bambu/done.txt",
        "results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam.bai",
        reads="results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam",
        gencode_transcriptome="results/setup/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/quantify_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.rds",
    log:
        "logs/quantify_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    benchmark:
        repeat(
            "benchmarks/quantify_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
            config["repetitions"],
        )
    conda:
        "../envs/r/bambu.yaml"
    threads: config["quantify_threads"]
    script:
        "../scripts/r/run_bambu_single_cell.R"


rule quantify_sc_run_isosceles:
    input:
        "results/setup/install_isosceles/done.txt",
        "results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam.bai",
        reads="results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam",
        gencode_transcriptome="results/setup/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/quantify_sc_run_isosceles/{tech}/{subsample_number}_{number_to_sample}/{sample}.rds",
    log:
        "logs/quantify_sc_run_isosceles/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    benchmark:
        repeat(
            "benchmarks/quantify_sc_run_isosceles/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
            config["repetitions"],
        )
    conda:
        "../envs/r/isosceles.yaml"
    threads: config["quantify_threads"]
    script:
        "../scripts/r/run_isosceles_sc.R"


rule quantify_kallisto_sc_lr:
    input:
        reads="results/preprocess/prepare_kallisto_sc_{tech}_fastq/{subsample_number}_{number_to_sample}/{sample}.fastq.gz",
        gencode_idx=f"results/setup/create_lr_kallisto_index/gencode_k-{config['lr_kallisto_index_k']}.idx",
        gencode_transcriptome_gmap_headered="results/setup/standardize_gtf_files/gencode_map_headered.txt",
    output:
        "results/quantify_kallisto_sc_lr/{tech}/{subsample_number}_{number_to_sample}/{sample}/matrix.abundance.mtx",
    log:
        "logs/quantify_kallisto_sc_lr/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    benchmark:
        "benchmarks/quantify_kallisto_sc_lr/{tech}/{subsample_number}_{number_to_sample}/{sample}.log"
    conda:
        "../envs/standalone/kallisto.yaml"
    threads: config["quantify_threads"]
    params:
        outdir="results/quantify_kallisto_sc_lr/{tech}/{subsample_number}_{number_to_sample}/{sample}",
        bus_threshold=config["kallisto_bus_threshold"],
        lr_kallisto_index_k=config["lr_kallisto_index_k"],
        compression_level=config["compression_level"],
        platform=lambda wildcards: (
            "PacBio"
            if wildcards.tech == "pb"
            else "ONT" if wildcards.tech == "ont" else ""
        ),
    shell:
        """
        conda list >> {log};

        kallisto bus \
            -t {threads} \
            --long \
            --threshold {params.bus_threshold} \
            -x "0,0,16:-1,0,0:0,16,0" \
            -i {input.gencode_idx} \
            -o {params.outdir} {input.reads} >> {log} 2>&1

        bustools sort -t {threads} {params.outdir}/output.bus \
            -o {params.outdir}/sorted.bus >> {log} 2>&1

        bustools count {params.outdir}/sorted.bus \
            -t {params.outdir}/transcripts.txt \
            -e {params.outdir}/matrix.ec \
            -g {input.gencode_transcriptome_gmap_headered} \
            -o {params.outdir}/count --cm -m >> {log} 2>&1

        kallisto quant-tcc \
            -t {threads} \
            --long -P {params.platform} \
            -f {params.outdir}/flens.txt \
            {params.outdir}/count.mtx \
            -i {input.gencode_idx} \
            -e {params.outdir}/count.ec.txt \
            -o {params.outdir} \
            >> {log} 2>&1
        """


rule quantify_sc_run_simpleaf:
    input:
        "results/setup/create_simpleaf_index/index/piscem_idx.sshash",
        "results/setup/create_simpleaf_index/index/piscem_idx.refinfo",
        first_reads="results/downsample_sc/illumina/{subsample_number}_{number_to_sample}/{sample}_R1.fastq.gz",
        second_reads="results/downsample_sc/illumina/{subsample_number}_{number_to_sample}/{sample}_R2.fastq.gz",
        t2g="results/setup/create_simpleaf_index/index/t2g_3col.tsv",
        index="results/setup/create_simpleaf_index",
    output:
        directory(
            "results/quantify_sc_run_simpleaf/{subsample_number}_{number_to_sample}/{sample}"
        ),
    log:
        "logs/quantify_sc_run_simpleaf/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/simpleaf.yaml"
    threads: config["quantify_threads"]
    params:
        index_dir="results/setup/create_simpleaf_index/index/piscem_idx",
    shell:
        """
        conda list &> {log};
        export ALEVIN_FRY_HOME="{input.index}";
        simpleaf set-paths &>> {log};
        ulimit -n 2048 &>> {log};
        simpleaf quant \
            --reads1 {input.first_reads} \
            --reads2 {input.second_reads} \
            --threads {threads} \
            --index {params.index_dir} \
            --t2g-map {input.t2g} \
            --chemistry 10xv3 --resolution cr-like \
            --unfiltered-pl --anndata-out \
            --output {output}
        """


rule quantify_sc_run_simpleaf_full_depth:
    input:
        "results/setup/create_simpleaf_index/index/piscem_idx.sshash",
        "results/setup/create_simpleaf_index/index/piscem_idx.refinfo",
        first_reads="/home/dwissel/data/roche_fmr1/illumina_sc/{sample}_R1.fastq.gz",
        second_reads="/home/dwissel/data/roche_fmr1/illumina_sc/{sample}_R2.fastq.gz",
        t2g="results/setup/create_simpleaf_index/index/t2g_3col.tsv",
        index="results/setup/create_simpleaf_index",
    output:
        directory("results/quantify_sc_run_simpleaf_full_depth/{sample}"),
    log:
        "logs/quantify_sc_run_simpleaf_full_depth/{sample}.log",
    conda:
        "../envs/standalone/simpleaf.yaml"
    threads: config["quantify_threads"]
    params:
        index_dir="results/setup/create_simpleaf_index/index/piscem_idx",
    shell:
        """
        conda list &> {log};
        export ALEVIN_FRY_HOME="{input.index}";
        simpleaf set-paths &>> {log};
        ulimit -n 2048 &>> {log};
        simpleaf quant \
            --reads1 {input.first_reads} \
            --reads2 {input.second_reads} \
            --threads {threads} \
            --index {params.index_dir} \
            --t2g-map {input.t2g} \
            --chemistry 10xv3 --resolution cr-like \
            --unfiltered-pl --anndata-out \
            --output {output}
        """


rule dump_isosceles_sc:
    input:
        rds="results/quantify_sc_run_isosceles/{tech}/{subsample_number}_{number_to_sample}/{sample}.rds",
        t2g="results/setup/standardize_gtf_files/gencode_map_headered.txt",
    output:
        mtx="results/quantify_sc_run_isosceles/{tech}/{subsample_number}_{number_to_sample}/{sample}/count.mtx",
        fts="results/quantify_sc_run_isosceles/{tech}/{subsample_number}_{number_to_sample}/{sample}/features.txt",
        bcs="results/quantify_sc_run_isosceles/{tech}/{subsample_number}_{number_to_sample}/{sample}/barcodes.txt",
    log:
        "results/quantify_sc_run_isosceles/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/r/formatting_env.yaml"
    threads: 6
    script:
        "../scripts/r/dump_isosceles_sc.R"


rule dump_bambu_sc:
    input:
        rds="results/quantify_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.rds",
    output:
        mtx="results/quantify_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}/count.mtx",
        fts="results/quantify_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}/features.txt",
        bcs="results/quantify_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}/barcodes.txt",
    log:
        "logs/quantify_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/r/formatting_env.yaml"
    threads: 6
    script:
        "../scripts/r/dump_bambu_sc.R"


SAMPLES = ["E3-1", "E3-2", "E3-3", "isoB11-1", "isoB11-2", "isoB11-3"]


def get_mtx_inputs(wildcards):
    if wildcards.method.lower() == "kallisto":
        return expand(
            "results/quantify_kallisto_sc_lr/{tech}/{subsample_number}_{number_to_sample}/{sample}/matrix.abundance.mtx",
            sample=SAMPLES,
            **wildcards,
        )

    elif wildcards.method.lower() == "oarfish":
        return expand(
            "results/quantify_sc_run_oarfish/{tech}/{subsample_number}_{number_to_sample}/{sample}.count.mtx",
            sample=SAMPLES,
            **wildcards,
        )

    else:
        return expand(
            "results/quantify_sc_run_{method}/{tech}/{subsample_number}_{number_to_sample}/{sample}/count.mtx",
            sample=SAMPLES,
            **wildcards,
        )


rule format_sce_lr:
    input:
        mtx=get_mtx_inputs,
        t2g="results/setup/standardize_gtf_files/gencode_map_headered.txt",
    output:
        tx_sce="results/format_quantify_subsampled/gencode/{tech}/{subsample_number}_{number_to_sample}/{method}/sce/tx_sce.rds",
        gene_sce="results/format_quantify_subsampled/gencode/{tech}/{subsample_number}_{number_to_sample}/{method}/sce/gene_sce.rds",
        tx_pb="results/format_quantify_subsampled/gencode/{tech}/{subsample_number}_{number_to_sample}/{method}/pseudobulk/transcript_counts_formatted.tsv",
        gene_pb="results/format_quantify_subsampled/gencode/{tech}/{subsample_number}_{number_to_sample}/{method}/pseudobulk/gene_counts_formatted.tsv",
    log:
        "logs/format_quantify_subsampled/gencode/{tech}/{subsample_number}_{number_to_sample}/{method}.log",
    conda:
        "../envs/r/formatting_env.yaml"
    threads: 4
    script:
        "../scripts/r/format_sce_lr.R"
