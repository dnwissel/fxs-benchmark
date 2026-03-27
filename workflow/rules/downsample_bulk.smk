configfile: "config/config.yaml"


rule downsample_short_read_fastq:
    input:
        first_reads="results/prepare_bulk/illumina_convert_mapped_bams_to_fastq/{sample}/{data_type}/{sample}.reads_1.fastq.gz",
        second_reads="results/prepare_bulk/illumina_convert_mapped_bams_to_fastq/{sample}/{data_type}/{sample}.reads_2.fastq.gz",
    output:
        first_reads="results/downsample_bulk/short_read_fastq/{subsample_number}_{read_number}_{data_type}/{sample}-r1.fastq.gz",
        second_reads="results/downsample_bulk/short_read_fastq/{subsample_number}_{read_number}_{data_type}/{sample}-r2.fastq.gz",
    log:
        "logs/downsample_bulk/short_read_fastq/{data_type}/{subsample_number}/{read_number}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    threads: config["align_threads"]
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.read_number,
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


rule downsample_long_read_fastq:
    input:
        "results/prepare_bulk/ont_and_pb_convert_mapped_bams_to_fastq/{data_type}/{sample}/{type}/{sample}.fastq.gz",
    output:
        "results/downsample_bulk/long_read_fastq/{type}/{subsample_number}_{read_number}_{data_type}/{sample}-r1.fastq.gz",
    log:
        "logs/downsample_bulk/long_read_fastq/{type}/{subsample_number}_{read_number}_{data_type}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    threads: config["align_threads"]
    params:
        seed=lambda wc: wc.subsample_number,
        number_to_sample=lambda wc: wc.read_number,
    shell:
        """
        conda list &> {log};
        seqtk sample -2 -s {params.seed} {input} \
            {params.number_to_sample} 2>> {log} | gzip > \
            {output} 2>> {log}
        """
