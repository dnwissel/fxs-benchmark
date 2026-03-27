rule fix_ont_headers:
    input:
        f"{config['sc_ont_directory']}/{{sample}}.fastq.gz",
    output:
        "results/downsample_sc/ont/{subsample_number}_{number_to_sample}/{sample}.fastq.gz",
    log:
        "logs/downsample_sc/ont/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    threads: config["stall_io_threads"]
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.number_to_sample,
    shell:
        """
        conda list &> {log};
        seqtk sample -2 -s {params.seed} {input} \
            {params.number_to_sample} 2>> {log} | gzip > \
            {output} 2>> {log}
        """


rule downsample_sc_ont:
    input:
        f"{config['sc_ont_directory']}/{{sample}}.fastq.gz",
    output:
        "results/downsample_sc/ont/{subsample_number}_{number_to_sample}/{sample}.fastq.gz",
    log:
        "logs/downsample_sc/ont/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    threads: config["stall_io_threads"]
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.number_to_sample,
    shell:
        """
        conda list &> {log};
        seqtk sample -2 -s {params.seed} {input} \
            {params.number_to_sample} 2>> {log} | pigz -p {threads} > \
            {output} 2>> {log}
        """


rule downsample_sc_illumina:
    input:
        first_reads=f"{config['sc_ill_directory']}/{{sample}}_R1.fastq.gz",
        second_reads=f"{config['sc_ill_directory']}/{{sample}}_R2.fastq.gz",
    output:
        first_reads="results/downsample_sc/illumina/{subsample_number}_{number_to_sample}/{sample}_R1.fastq.gz",
        second_reads="results/downsample_sc/illumina/{subsample_number}_{number_to_sample}/{sample}_R2.fastq.gz",
    log:
        "logs/downsample_sc/illumina/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    threads: config["align_threads"]
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.number_to_sample,
    shell:
        """
        conda list &> {log};
        seqtk sample -2 -s {params.seed} {input.first_reads} \
            {params.number_to_sample} 2>> {log} | gzip > \
            {output.first_reads} 2>> {log};
        seqtk sample -2 -s {params.seed} {input.second_reads} \
            {params.number_to_sample} 2>> {log} | gzip > \
            {output.second_reads} 2>> {log}
        """


def get_merge_inputs(wildcards):
    inputs = [f"{config['sc_pb_directory']}/{wildcards.sample}_run1-segmented.bam"]

    if wildcards.sample != "isoB11-3":
        inputs.append(
            f"{config['sc_pb_directory']}/{wildcards.sample}_run2-segmented.bam"
        )

    return inputs


rule merge_pb_sc_runs:
    input:
        get_merge_inputs,
    output:
        "results/merge_pb_sc_runs/{sample}-segmented.bam",
    log:
        "logs/merge_pb_sc_runs/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["stall_io_threads"]
    shell:
        """
        conda list &> {log};
        samtools cat -o {output} {input} &>> {log}
        """


rule subsample_bam_fixed_number_reformat:
    input:
        "results/merge_pb_sc_runs/{sample}-segmented.bam",
    output:
        "results/downsample_sc/pb/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/downsample_sc/pb/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/bbtools.yaml"
    threads: config["stall_io_threads"]
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.number_to_sample,
    shell:
        """
        conda list &> {log};
        reformat.sh \
            in={input} \
            out={output} \
            samplereadstarget={params.number_to_sample} \
            sampleseed={params.seed} \
            2>> {log}
        """
