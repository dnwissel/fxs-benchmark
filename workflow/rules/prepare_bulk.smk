
configfile: "config/config.yaml"


rule illumina_run_trim_galore:
    input:
        read_1=f"{config['bulk_ill_directory']}/{{sample}}_R1.fastq.gz",
        read_2=f"{config['bulk_ill_directory']}/{{sample}}_R2.fastq.gz",
    output:
        reads_1="results/prepare_bulk/illumina_run_trim_galore/{sample}/{sample}_R1_val_1.fq.gz",
        reads_2="results/prepare_bulk/illumina_run_trim_galore/{sample}/{sample}_R2_val_2.fq.gz",
    log:
        "logs/prepare_bulk/illumina_run_trim_galore/{sample}.out",
    conda:
        "../envs/standalone/trim_galore.yaml"
    params:
        min_quality=config["bulk_ill_trim_galore_min_q"],
        min_length=config["bulk_ill_trim_galore_min_length"],
        output_dir="results/prepare_bulk/illumina_run_trim_galore/{sample}",
    shell:
        """
        conda list &> {log};
        trim_galore -q {params.min_quality} --phred33 \
            --length {params.min_length} -o {params.output_dir} \
            --paired {input.read_1} {input.read_2}  &>> {log}
        """


rule illumina_run_star_genome:
    input:
        reads_1="results/prepare_bulk/illumina_run_trim_galore/{sample}/{sample}_R1_val_1.fq.gz",
        reads_2="results/prepare_bulk/illumina_run_trim_galore/{sample}/{sample}_R2_val_2.fq.gz",
        genome="results/setup/index_star",
    output:
        "results/prepare_bulk/illumina_run_star_genome/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    log:
        "logs/prepare_bulk/run_star_genome/{sample}.out",
    conda:
        "../envs/standalone/star.yaml"
    threads: config["align_threads"]
    params:
        "results/prepare_bulk/illumina_run_star_genome/{sample}/{sample}_",
    shell:
        """
        conda list &> {log};
        STAR --genomeDir {input.genome} --readFilesIn {input.reads_1} \
            {input.reads_2} --runThreadN {threads} \
            --outFileNamePrefix {params} \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesCommand gunzip -c \
            --outSAMstrandField intronMotif &>> {log}
        """


rule illumina_index_star_genome:
    input:
        "results/prepare_bulk/illumina_run_star_genome/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "results/prepare_bulk/illumina_run_star_genome/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
    log:
        "logs/prepare_bulk/index_illumina_star/{sample}.out",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        conda list &> {log};
        sleep 61s;
        samtools index {input} &>> {log}
        """


rule illumina_separate_genome_mappings:
    input:
        bam="results/prepare_bulk/illumina_run_star_genome/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        index="results/prepare_bulk/illumina_run_star_genome/{sample}/{sample}_Aligned.sortedByCoord.out.bam.bai",
    output:
        gencode_reads="results/prepare_bulk/illumina_separate_genome_mappings/{sample}/{sample}.aligned.gencode.sorted.bam",
        gencode_ix="results/prepare_bulk/illumina_separate_genome_mappings/{sample}/{sample}.aligned.gencode.sorted.bam.bai",
        sirv_reads="results/prepare_bulk/illumina_separate_genome_mappings/{sample}/{sample}.aligned.sirv.sorted.bam",
        sirv_ix="results/prepare_bulk/illumina_separate_genome_mappings/{sample}/{sample}.aligned.sirv.sorted.bam.bai",
        ercc_reads="results/prepare_bulk/illumina_separate_genome_mappings/{sample}/{sample}.aligned.ercc.sorted.bam",
        ercc_ix="results/prepare_bulk/illumina_separate_genome_mappings/{sample}/{sample}.aligned.ercc.sorted.bam.bai",
    log:
        "logs/prepare_bulk/illumina_separate_genome_mappings/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["stall_io_threads"]
    params:
        lines_to_cut=1,
        gencode_chromosomes="chr",
        sirv_chromosomes="SIRV",
        ercc_chromosomes="ERCC",
    shell:
        """
        conda list &> {log};
        samtools idxstats {input.bam} | \
            cut -f {params.lines_to_cut} | \
            grep {params.gencode_chromosomes}  | \
            xargs samtools view -b {input.bam} > {output.gencode_reads};
        samtools idxstats {input.bam} | \
            cut -f {params.lines_to_cut} | \
            grep {params.sirv_chromosomes} | \
            xargs samtools view -F 0x904 -b {input.bam} \
            > {output.sirv_reads};
        samtools idxstats {input.bam} | \
            cut -f {params.lines_to_cut} | \
            grep {params.ercc_chromosomes} | \
            xargs samtools view -F 0x904 -b {input.bam} \
            > {output.ercc_reads};
        sleep 61s &>> {log};
        samtools index {output.gencode_reads} &>> {log};
        samtools index {output.sirv_reads} &>> {log};
        samtools index {output.ercc_reads} &>> {log}
        """


rule illumina_convert_mapped_bams_to_fastq:
    input:
        "results/prepare_bulk/illumina_separate_genome_mappings/{sample}/{sample}.aligned.{type}.sorted.bam",
    output:
        singletons="results/prepare_bulk/illumina_convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.singletons.fastq.gz",
        first_reads="results/prepare_bulk/illumina_convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_1.fastq.gz",
        second_reads="results/prepare_bulk/illumina_convert_mapped_bams_to_fastq/{sample}/{type}/{sample}.reads_2.fastq.gz",
    log:
        "logs/prepare_bulk/illumina_convert_mapped_bams_to_fastq/{type}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["stall_io_threads"]
    params:
        memory=config["sort_bam_memory_gb"],
        compression_level=config["compression_level"],
    shell:
        """
        conda list &> {log};
        samtools sort -n -m{params.memory}g -@ {threads} {input} \
            2>> {log} | samtools fastq -@ {threads} \
            -1 {output.first_reads} -2 {output.second_reads} \
            -s {output.singletons} -c {params.compression_level} \
            - &>> {log}
        """


rule ont_run_pychopper:
    input:
        read=f"{config['bulk_ont_directory']}/{{sample}}_R1.fastq.gz",
    output:
        unclassified_reads="results/prepare_bulk/ont_run_pychopper/{sample}/unclassified.fq.gz",
        rescued_reads="results/prepare_bulk/ont_run_pychopper/{sample}/rescued.fq.gz",
        full_length_reads="results/prepare_bulk/ont_run_pychopper/{sample}/full_length_output.fq.gz",
    log:
        "logs/preprocessing/pychopper/{sample}/log.out",
    conda:
        "../envs/standalone/pychopper.yaml"
    threads: config["stall_io_threads"]
    params:
        outdir="results/prepare_bulk/ont_run_pychopper/{sample}",
        kit=config["bulk_ont_pychopper_kit"],
        y=config["bulk_ont_pychopper_y"],
    shell:
        """
        conda list &> {log};
        pychopper -Y {params.y} -k {params.kit} \
            -r {params.outdir}/report.pdf -t {threads} \
            -u {params.outdir}/unclassified.fq \
            -w {params.outdir}/rescued.fq {input.read} \
            {params.outdir}/full_length_output.fq &>> {log};
        gzip {params.outdir}/unclassified.fq;
        gzip {params.outdir}/rescued.fq;
        gzip {params.outdir}/full_length_output.fq
        """


rule ont_run_minimap2_genome:
    input:
        read="results/prepare_bulk/ont_run_pychopper/{sample}/full_length_output.fq.gz",
        transcriptome="results/setup/convert_gtfs_to_beds/overall.bed",
        genome="results/setup/concatenate_genomes/genome.fa",
    output:
        "results/prepare_bulk/ont_run_minimap2_genome/{sample}/{sample}.aligned.bam",
    log:
        "logs/prepare_bulk/ont_run_minimap2_genome/{sample}.stdout",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: config["map_threads"]
    params:
        map_threads=config["align_threads"],
        sort_threads=config["sort_threads"],
        sort_memory_gb=config["sort_memory_gb"],
    shell:
        """
        conda list &> {log};
        minimap2 -ax splice --junc-bed {input.transcriptome} \
            -t {params.map_threads} \
            {input.genome} {input.read} | samtools sort \
            -@ {params.sort_threads} \
            -m{params.sort_memory_gb}g \
            -o {output}##idx##{output}.bai --write-index - &>> {log}
        """


rule pb_convert_to_bam:
    input:
        read=f"{config['bulk_pb_directory']}/{{sample}}.bam",
    output:
        "results/prepare_bulk/pb_convert_to_bam/{sample}.fastq.gz",
    log:
        "logs/prepare_bulk/pb_convert_to_bam/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: 6
    shell:
        """
        conda list &> {log};
        samtools bam2fq -@ {threads} {input.read} | gzip > \
            {output} 2>> {log}
        """


rule pb_run_minimap2_genome:
    input:
        read="results/prepare_bulk/pb_convert_to_bam/{sample}.fastq.gz",
        transcriptome="results/setup/convert_gtfs_to_beds/overall.bed",
        genome="results/setup/concatenate_genomes/genome.fa",
    output:
        "results/prepare_bulk/pb_run_minimap2_genome/{sample}/{sample}.aligned.bam",
    log:
        "logs/prepare_bulk/pb_run_minimap2_genome/{sample}.stdout",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: config["map_threads"]
    params:
        map_threads=config["align_threads"],
        sort_threads=config["sort_threads"],
        sort_memory_gb=config["sort_memory_gb"],
    shell:
        """
        conda list &> {log};
        minimap2 -ax splice:hq --junc-bed {input.transcriptome} \
            -t {params.map_threads} \
            {input.genome} {input.read} | samtools sort \
            -@ {params.sort_threads} \
            -m{params.sort_memory_gb}g \
            -o {output}##idx##{output}.bai --write-index - &>> {log}
        """


rule ont_and_pb_separate_genome_mappings:
    input:
        bam="results/prepare_bulk/{type}_run_minimap2_genome/{sample}/{sample}.aligned.bam",
    output:
        gencode_reads="results/prepare_bulk/{type}_separate_genome_mappings/{sample}/{sample}.aligned.gencode.sorted.bam",
        gencode_ix="results/prepare_bulk/{type}_separate_genome_mappings/{sample}/{sample}.aligned.gencode.sorted.bam.bai",
        sirv_reads="results/prepare_bulk/{type}_separate_genome_mappings/{sample}/{sample}.aligned.sirv.sorted.bam",
        sirv_ix="results/prepare_bulk/{type}_separate_genome_mappings/{sample}/{sample}.aligned.sirv.sorted.bam.bai",
        ercc_reads="results/prepare_bulk/{type}_separate_genome_mappings/{sample}/{sample}.aligned.ercc.sorted.bam",
        ercc_ix="results/prepare_bulk/{type}_separate_genome_mappings/{sample}/{sample}.aligned.ercc.sorted.bam.bai",
    log:
        "logs/prepare_bulk/{type}_ont_and_pb_separate_genome_mappings/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["stall_io_threads"]
    params:
        lines_to_cut=1,
        gencode_chromosomes="chr",
        sirv_chromosomes="SIRV",
        ercc_chromosomes="ERCC",
    shell:
        """
        conda list &> {log};
        samtools idxstats {input.bam} | \
            cut -f {params.lines_to_cut} | \
            grep {params.gencode_chromosomes}  | \
            xargs samtools view -b {input.bam} > {output.gencode_reads};
        samtools idxstats {input.bam} | \
            cut -f {params.lines_to_cut} | \
            grep {params.sirv_chromosomes} | \
            xargs samtools view -F 0x904 -b {input.bam} \
            > {output.sirv_reads};
        samtools idxstats {input.bam} | \
            cut -f {params.lines_to_cut} | \
            grep {params.ercc_chromosomes} | \
            xargs samtools view -F 0x904 -b {input.bam} \
            > {output.ercc_reads};
        sleep 61s &>> {log};
        samtools index {output.gencode_reads} &>> {log};
        samtools index {output.sirv_reads} &>> {log};
        samtools index {output.ercc_reads} &>> {log}
        """


rule ont_and_pb_convert_mapped_bams_to_fastq:
    input:
        "results/prepare_bulk/{type}_separate_genome_mappings/{sample}/{sample}.aligned.{data_type}.sorted.bam",
    output:
        "results/prepare_bulk/ont_and_pb_convert_mapped_bams_to_fastq/{data_type}/{sample}/{type}/{sample}.fastq.gz",
    log:
        "logs/prepare_bulk/ont_and_pb_convert_mapped_bams_to_fastq/{type}/{data_type}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["align_threads"]
    params:
        memory=config["sort_bam_memory_gb"],
        compression_level=config["compression_level"],
    shell:
        """
        conda list &> {log};
        samtools sort -n -m{params.memory}g -@ 4 {input} \
            2>> {log} | samtools fastq -@ 4 \
            -c {params.compression_level} - > {output} 2>> {log}
        """
