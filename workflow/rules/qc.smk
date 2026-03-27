configfile: "config/config.yaml"


rule star_solo_short_read:
    input:
        reads_1="/home/dwissel/data/roche_fmr1/illumina_sc/{sample}_R1.fastq.gz",
        reads_2="/home/dwissel/data/roche_fmr1/illumina_sc/{sample}_R2.fastq.gz",
        star_index="results/setup/index_star",
        gtf="results/setup/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
        whitelist="config/3M-february-2018.txt",
    output:
        bam="results/qc/star_solo/{sample}/Aligned.sortedByCoord.out.bam",
    log:
        "logs/qc/star_solo/{sample}.log",
    conda:
        "../envs/standalone/star_and_samtools.yaml"
    threads: config["align_threads"]
    params:
        out_prefix="results/qc/star_solo/{sample}/",
        cb_start=1,
        cb_len=16,
        umi_start=17,
        umi_len=12,
    shell:
        """
        mkdir -p {params.out_prefix}

        STAR --runThreadN {threads} \
             --genomeDir {input.star_index} \
             --readFilesIn {input.reads_2} {input.reads_1} \
             --readFilesCommand zcat \
             --sjdbGTFfile {input.gtf} \
             --outFileNamePrefix {params.out_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
             --soloType CB_UMI_Simple \
             --soloCBwhitelist {input.whitelist} \
             --soloCBstart {params.cb_start} --soloCBlen {params.cb_len} \
             --soloUMIstart {params.umi_start} --soloUMIlen {params.umi_len} \
             --soloUMIdedup 1MM_All \
             --soloUMIfiltering - \
             --outSAMunmapped Within \
             &> {log}
        """


rule downsample_long_read_fastq_qc_ont:
    input:
        lambda wildcards: (
            f"{config['bulk_ont_directory']}/{wildcards.sample}_R1.fastq.gz"
            if wildcards.platform == "bulk"
            else f"{config['sc_ont_directory']}/{wildcards.sample}.fastq.gz"
        ),
    output:
        "results/qc/downsample_long_read_fastq_qc_ont/{read_number}/{platform}/{sample}.fastq.gz",
    log:
        "logs/qc/downsample_long_read_fastq_qc_ont/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    threads: 6
    params:
        seed=config["seed"],
        number_to_sample=lambda wc: wc.read_number,
    shell:
        """
        conda list &> {log};
        seqtk sample -2 -s {params.seed} {input} \
            {params.number_to_sample} 2>> {log} | gzip > \
            {output} 2>> {log}
        """


rule downsample_long_read_fastq_qc_pb:
    input:
        lambda wildcards: (
            f"{config['bulk_pb_directory']}/{wildcards.sample}.bam"
            if wildcards.platform == "bulk"
            else "results/merge_pb_sc_runs/{sample}-segmented.bam"
        ),
    output:
        "results/qc/downsample_long_read_fastq_qc_pb/{read_number}/{platform}/{sample}.fastq.gz",
    log:
        "logs/qc/downsample_long_read_fastq_qc_pb/{read_number}/{platform}/{sample}.log",
    wildcard_constraints:
        platform="bulk",
    conda:
        "../envs/standalone/bbtools.yaml"
    threads: 6
    params:
        seed=config["seed"],
        number_to_sample=lambda wc: wc.read_number,
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


rule downsample_long_read_fastq_qc_pb_sc:
    input:
        "results/merge_pb_sc_runs/{sample}-segmented.bam",
    output:
        "results/qc/downsample_long_read_fastq_qc_pb_sc/{read_number}/{platform}/{sample}.bam",
    log:
        "logs/qc/downsample_long_read_fastq_qc_pb_sc/{read_number}/{platform}/{sample}.log",
    wildcard_constraints:
        platform="single_cell",
    conda:
        "../envs/standalone/bbtools.yaml"
    threads: 6
    params:
        seed=config["seed"],
        number_to_sample=lambda wc: wc.read_number,
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


rule preprocess_long_read_fastq_qc_pb_sc:
    input:
        reads="results/qc/downsample_long_read_fastq_qc_pb_sc/{read_number}/{platform}/{sample}.bam",
        primers="config/pb_3p_primers.fasta",
        barcodes="config/3M-february-2018-REVERSE-COMPLEMENTED.txt",
    output:
        fastq="results/qc/downsample_long_read_fastq_qc_pb/{read_number}/{platform}/{sample}.fastq.gz",
    log:
        "logs/qc/preprocess_long_read_fastq_qc_pb_sc/{read_number}/{platform}/{sample}.log",
    wildcard_constraints:
        platform="single_cell",
    conda:
        "../envs/standalone/isoseq_bam.yaml"
    threads: 6
    params:
        dir="results/qc/downsample_long_read_fastq_qc_pb/{read_number}/{platform}",
        tmp_bam="results/qc/downsample_long_read_fastq_qc_pb/{read_number}/{platform}/{sample}.fltnc.tmp.bam",
        out_prefix="results/qc/downsample_long_read_fastq_qc_pb/{read_number}/{platform}/{sample}",
        outfile_primers="results/qc/downsample_sc/preprocess_pb/{read_number}/{platform}/{sample}.bam",
        primers="results/qc/downsample_sc/preprocess_pb/{read_number}/{platform}/{sample}.5p--3p.bam",
        tag="results/qc/downsample_sc/preprocess_pb/{read_number}/{platform}/{sample}.fltn.bam",
    shell:
        """
        exec &> {log};
        
        mkdir -p {params.dir};
        mkdir -p $(dirname {params.outfile_primers}); 

        lima {input.reads} {input.primers} {params.outfile_primers} --no-reports \
            --num-threads {threads} --isoseq;

        isoseq tag {params.primers} {params.tag} \
            --design T-12U-16B --num-threads {threads};

        rm {params.primers};
        isoseq refine {params.tag} {input.primers} {params.tmp_bam} \
            --require-polya --num-threads {threads};

        bam2fastq -o {params.out_prefix} \
            -j {threads} \
            -c 6 \
            {params.tmp_bam};

        rm {params.tmp_bam} {params.tag};
        """


rule downsample_short_read_fastq_qc:
    input:
        first_reads=lambda wildcards: (
            f"{config['bulk_ill_directory']}/{wildcards.sample}_R1.fastq.gz"
            if wildcards.platform == "bulk"
            else f"{config['sc_ill_directory']}/{wildcards.sample}_R1.fastq.gz"
        ),
        second_reads=lambda wildcards: (
            f"{config['bulk_ill_directory']}/{wildcards.sample}_R2.fastq.gz"
            if wildcards.platform == "bulk"
            else f"{config['sc_ill_directory']}/{wildcards.sample}_R2.fastq.gz"
        ),
    output:
        first_reads="results/qc/downsample_short_read_fastq_qc_ill/{read_number}/{platform}/{sample}_R1.fastq.gz",
        second_reads="results/qc/downsample_short_read_fastq_qc_ill/{read_number}/{platform}/{sample}_R2.fastq.gz",
    log:
        "logs/qc/downsample_short_read_fastq_qc_ill/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    threads: 6
    params:
        seed=config["seed"],
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


rule qc_calculate_read_lengths_tsv_short_reads:
    input:
        first_reads="results/qc/downsample_short_read_fastq_qc_ill/{read_number}/{platform}/{sample}_R1.fastq.gz",
        second_reads="results/qc/downsample_short_read_fastq_qc_ill/{read_number}/{platform}/{sample}_R2.fastq.gz",
    output:
        first="results/qc/qc_calculate_read_lengths_tsv_short_reads/{read_number}/{platform}/{sample}_R1.txt",
        second="results/qc/qc_calculate_read_lengths_tsv_short_reads/{read_number}/{platform}/{sample}_R2.txt",
    log:
        "logs/qc/qc_calculate_read_lengths_tsv_short_reads/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{print $name"\\t"length($seq)}}' \
            {input.first_reads} > {output.first} 2>> {log};
        bioawk -c fastx '{{print $name"\\t"length($seq)}}' \
            {input.second_reads} > {output.second} 2>> {log}
        """


rule qc_calculate_base_quality_tsv_short_reads:
    input:
        first_reads="results/qc/downsample_short_read_fastq_qc_ill/{read_number}/{platform}/{sample}_R1.fastq.gz",
        second_reads="results/qc/downsample_short_read_fastq_qc_ill/{read_number}/{platform}/{sample}_R2.fastq.gz",
    output:
        first="results/qc/qc_calculate_base_quality/{read_number}/{platform}/{sample}_R1.txt",
        second="results/qc/qc_calculate_base_quality/{read_number}/{platform}/{sample}_R2.txt",
    log:
        "logs/qc/qc_calculate_base_quality_tsv_short_reads/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{print $name"\\t"meanqual($qual)}}' \
            {input.first_reads} > {output.first} 2>> {log};
        bioawk -c fastx '{{print $name"\\t"meanqual($qual)}}' \
            {input.second_reads} > {output.second} 2>> {log};
        """


rule qc_calculate_read_lengths_tsv:
    input:
        "results/qc/downsample_long_read_fastq_qc_{tech}/{read_number}/{platform}/{sample}.fastq.gz",
    output:
        "results/qc/qc_calculate_read_lengths/{tech}/{read_number}/{platform}/{sample}.txt",
    log:
        "logs/qc/qc_calculate_read_lengths/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{print $name"\\t"length($seq)}}' \
            {input} > {output} 2>> {log}
        """


rule qc_calculate_base_quality_tsv:
    input:
        "results/qc/downsample_long_read_fastq_qc_{tech}/{read_number}/{platform}/{sample}.fastq.gz",
    output:
        "results/qc/qc_calculate_base_quality/{tech}/{read_number}/{platform}/{sample}.txt",
    log:
        "logs/qc/qc_calculate_base_quality/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/qc.yaml"
    shell:
        """
        conda list > {log};
        bioawk -c fastx '{{print $name"\\t"meanqual($qual)}}' \
            {input} > {output} 2>> {log}
        """


def get_fastq_for_star(wildcards):

    base_path = f"results/qc/downsample_short_read_fastq_qc_ill/{wildcards.read_number}/{wildcards.platform}/{wildcards.sample}"

    r1 = f"{base_path}_R1.fastq.gz"
    r2 = f"{base_path}_R2.fastq.gz"

    if wildcards.platform == "bulk":
        return [r1, r2]
    elif wildcards.platform == "sc":
        return [r2]
    else:
        return [r1, r2]


rule qc_align_genome_short_read:
    input:
        reads=get_fastq_for_star,
        star_index="results/setup/index_star",
        gtf="results/setup/standardize_gtf_files/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        bam="results/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}.bam",
        idx="results/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}.bam.bai",
        tx_bam="results/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}.transcriptome.bam",
        name_sorted_bam="results/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}.namesorted.bam",
    log:
        "logs/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/star_and_samtools.yaml"
    threads: config["align_threads"]
    params:
        out_prefix="results/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}_",
        memory_gb=config.get("align_sort_bam_memory_gb", 30),
    shell:
        """
        STAR --runThreadN {threads} \
             --genomeDir {input.star_index} \
             --readFilesIn {input.reads} \
             --readFilesCommand zcat \
             --outSAMattributes Standard NM \
             --sjdbGTFfile {input.gtf} \
             --outFileNamePrefix {params.out_prefix} \
             --outSAMtype BAM SortedByCoordinate \
             --quantMode TranscriptomeSAM \
             --outSAMunmapped Within \
             &> {log}

        mv {params.out_prefix}Aligned.sortedByCoord.out.bam {output.bam} >> {log} 2>&1
        samtools index {output.bam} >> {log} 2>&1

        samtools sort -n -@ {threads} -o {output.name_sorted_bam} {output.bam} >> {log} 2>&1

        mv {params.out_prefix}Aligned.toTranscriptome.out.bam {output.tx_bam} >> {log} 2>&1

        rm {params.out_prefix}Log.out {params.out_prefix}Log.progress.out {params.out_prefix}SJ.out.tab
        """


rule qc_align_genome:
    input:
        reads="results/qc/downsample_long_read_fastq_qc_{tech}/{read_number}/{platform}/{sample}.fastq.gz",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_bed="results/setup/convert_gtfs_to_beds/gencode.bed",
    output:
        "results/qc/qc_align_genome/{tech}/{read_number}/{platform}/{sample}.bam",
    log:
        "logs/qc/qc_align_genome/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: config["align_threads"]
    params:
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
    shell:
        """
        conda list &> {log};
        if [[ "{wildcards.tech}" == "pb" ]]; then
            minimap2 -ax splice:hq --junc-bed {input.gencode_bed} \
                -t {threads} {input.gencode_genome} \
                {input.reads} 2>> {log} | \
                samtools sort -@ {params.align_sort_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g -o {output} \
                - &>> {log};
            sleep 61s &>> {log};
            samtools index {output} &>> {log}
        elif [[ "{wildcards.tech}" == "ont" ]]; then
            minimap2 -ax splice --junc-bed {input.gencode_bed} \
                -t {threads} {input.gencode_genome} \
                {input.reads} 2>> {log} | \
                samtools sort -@ {params.align_sort_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g -o {output} \
                - &>> {log};
            sleep 61s &>> {log};
            samtools index {output} &>> {log}
        fi
        """


rule qc_align_transcriptome:
    input:
        reads="results/qc/downsample_long_read_fastq_qc_{tech}/{read_number}/{platform}/{sample}.fastq.gz",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
    output:
        "results/qc/qc_align_transcriptome/{tech}/{read_number}/{platform}/{sample}.bam",
    log:
        "logs/qc/qc_align_transcriptome/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: config["align_threads"]
    params:
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
        n_secondary_alignments=100,
    shell:
        """
        if [[ "{wildcards.tech}" == "pb" ]]; then
            preset="map-hifi"
        elif [[ "{wildcards.tech}" == "ont" ]]; then
            preset="lr:hq"
        fi

        minimap2 --eqx -N {params.n_secondary_alignments} -ax $preset \
            -t {threads} {input.gencode_transcriptome} \
            {input.reads} 2> {log} | \
            samtools view -@ {threads} -b -o {output} - 2>> {log}
        """


rule qc_add_alignment_md:
    input:
        reads="results/qc/qc_align_genome/{tech}/{read_number}/{platform}/{sample}.bam",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/qc/qc_add_alignment_md/{tech}/{read_number}/{platform}/{sample}.bam",
    log:
        "logs/qc/qc_add_alignment_md/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        conda list &> {log};
        samtools calmd -@ {threads} -b {input.reads} {input.gencode_genome} \
            > {output} 2>> {log}
        """


rule qc_add_alignment_md_short_reads:
    input:
        reads="results/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}.namesorted.bam",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/qc/qc_add_alignment_md_short_reads/{read_number}/{platform}/{sample}.bam",
    log:
        "logs/qc/qc_add_alignment_md_short_reads/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        conda list &> {log};
        samtools calmd -@ {threads} -b {input.reads} {input.gencode_genome} \
            > {output} 2>> {log}
        """


rule qc_calculate_alignment_qc:
    input:
        "results/qc/qc_add_alignment_md/{tech}/{read_number}/{platform}/{sample}.bam",
    output:
        "results/qc/qc_calculate_alignment_qc/{tech}/{read_number}/{platform}/{sample}.tsv",
    log:
        "logs/qc/qc_calculate_alignment_qc/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/py/pysam.yaml"
    params:
        junction_regex="[2-9][0-9][0-9]*N",
    script:
        "../scripts/py/calculate_alignment_qc_kinnex_rev.py"


rule qc_calculate_alignment_qc_short_read_single_cell:
    input:
        "results/qc/qc_add_alignment_md_short_reads/{read_number}/single_cell/{sample}.bam",
    output:
        "results/qc/qc_calculate_alignment_qc_short_read_single_cell/{read_number}/single_cell/{sample}.tsv",
    log:
        "logs/qc/qc_calculate_alignment_qc_short_read_single_cell/{read_number}/single_cell/{sample}.log",
    conda:
        "../envs/py/pysam.yaml"
    params:
        junction_regex="[2-9][0-9][0-9]*N",
    script:
        "../scripts/py/calculate_alignment_qc_kinnex_rev.py"


rule qc_calculate_alignment_qc_short_read_bulk:
    input:
        "results/qc/qc_add_alignment_md_short_reads/{read_number}/bulk/{sample}.bam",
    output:
        "results/qc/qc_calculate_alignment_qc_short_read_bulk/{read_number}/bulk/{sample}.tsv",
    log:
        "logs/qc/qc_calculate_alignment_qc_short_read_bulk/{read_number}/bulk/{sample}.log",
    conda:
        "../envs/py/pysam.yaml"
    params:
        junction_regex="[2-9][0-9][0-9]*N",
    script:
        "../scripts/py/calculate_alignment_qc_ill_rev.py"


rule qc_prepare_sampled_read_quality_frame:
    input:
        lengths="results/qc/qc_calculate_read_lengths/{tech}/{read_number}/{platform}/{sample}.txt",
        quality="results/qc/qc_calculate_base_quality/{tech}/{read_number}/{platform}/{sample}.txt",
        alignment="results/qc/qc_calculate_alignment_qc/{tech}/{read_number}/{platform}/{sample}.tsv",
    output:
        "results/qc/qc_prepare_sampled_read_quality_frame/{tech}/{read_number}/{platform}/{sample}.tsv",
    log:
        "logs/qc/qc_prepare_sampled_read_quality_frame/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/r/qual_qc.yaml"
    script:
        "../scripts/r/qc_prepare_sampled_read_quality_frame_kinnex_rev.R"


rule qc_prepare_sampled_read_quality_frame_short_read_single_cell:
    input:
        lengths="results/qc/qc_calculate_read_lengths_tsv_short_reads/{read_number}/single_cell/{sample}_R2.txt",
        quality="results/qc/qc_calculate_base_quality/{read_number}/single_cell/{sample}_R2.txt",
        alignment="results/qc/qc_calculate_alignment_qc_short_read_single_cell/{read_number}/single_cell/{sample}.tsv",
    output:
        "results/qc/qc_prepare_sampled_read_quality_frame_short_read_single_cell/{read_number}/single_cell/{sample}.tsv",
    log:
        "logs/qc/qc_prepare_sampled_read_quality_frame_short_read_single_cell/{read_number}/single_cell/{sample}.log",
    conda:
        "../envs/r/qual_qc.yaml"
    script:
        "../scripts/r/qc_prepare_sampled_read_quality_frame_kinnex_rev.R"


rule qc_prepare_sampled_read_quality_frame_short_read_bulk:
    input:
        first_lengths="results/qc/qc_calculate_read_lengths_tsv_short_reads/{read_number}/bulk/{sample}_R1.txt",
        second_lengths="results/qc/qc_calculate_read_lengths_tsv_short_reads/{read_number}/bulk/{sample}_R2.txt",
        first_quality="results/qc/qc_calculate_base_quality/{read_number}/bulk/{sample}_R1.txt",
        second_quality="results/qc/qc_calculate_base_quality/{read_number}/bulk/{sample}_R2.txt",
        alignment="results/qc/qc_calculate_alignment_qc_short_read_bulk/{read_number}/bulk/{sample}.tsv",
    output:
        "results/qc/qc_prepare_sampled_read_quality_frame_short_read_bulk/{read_number}/bulk/{sample}.tsv",
    log:
        "logs/qc/qc_prepare_sampled_read_quality_frame_short_read_bulk/{read_number}/bulk/{sample}.log",
    conda:
        "../envs/r/qual_qc.yaml"
    script:
        "../scripts/r/qc_prepare_sampled_read_quality_frame_illumina_rev.R"


rule qc_subsample_transcripts_three_prime_bias_binned:
    input:
        transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
    output:
        less_one_kb="results/qc/qc_subsample_transcripts_three_prime_bias_binned/lt1kb.gtf",
        one_two_kb="results/qc/qc_subsample_transcripts_three_prime_bias_binned/1-2kb.gtf",
        two_three_kb="results/qc/qc_subsample_transcripts_three_prime_bias_binned/2-3kb.gtf",
        three_six_kb="results/qc/qc_subsample_transcripts_three_prime_bias_binned/3-6kb.gtf",
        greater_six_kb="results/qc/qc_subsample_transcripts_three_prime_bias_binned/gt6kb.gtf",
    log:
        "logs/qc/qc_subsample_transcripts_three_prime_bias_binned/out.log",
    conda:
        "../envs/r/base.yaml"
    params:
        n_subsample=config["three_prime_bias_subsample"],
        seed=config["seed"],
    script:
        "../scripts/r/subsample_transcripts_three_prime_bias_binned.R"


rule qc_convert_gtf_to_bed_binned:
    input:
        "results/qc/qc_subsample_transcripts_three_prime_bias_binned/{nam}.gtf",
    output:
        "results/qc/convert_gtf_to_bed_binned/{nam}.bed",
    log:
        "logs/qc/convert_gtf_to_bed_binned/{nam}.log",
    conda:
        "../envs/standalone/ucsc_genepred.yaml"
    shell:
        """
        conda list &> "{log}";
        cat "{input}" 2>> "{log}" |\
            gtfToGenePred /dev/stdin /dev/stdout 2>> "{log}" |\
            genePredToBed /dev/stdin /dev/stdout > \
            "{output}" 2>> "{log}"
        """


rule qc_calculate_three_prime_bias_binned:
    input:
        gencode_bed="results/qc/convert_gtf_to_bed_binned/{nam}.bed",
        bam="results/qc/qc_align_genome/{tech}/{read_number}/{platform}/{sample}.bam",
    output:
        "results/qc/qc_calculate_three_prime_bias_binned/{tech}/{read_number}/{platform}/{sample}/{nam}/three_prime_bias.geneBodyCoverage.txt",
    log:
        "logs/qc/qc_calculate_three_prime_bias_binned/{tech}/{read_number}/{platform}/{sample}/{nam}.log",
    conda:
        "../envs/standalone/rseqc.yaml"
    threads: 1
    params:
        output_path="results/qc/qc_calculate_three_prime_bias_binned/{tech}/{read_number}/{platform}/{sample}/{nam}/three_prime_bias",
    shell:
        """
        conda list > "{log}";
        geneBody_coverage.py \
            -r "{input.gencode_bed}" \
            -o "{params.output_path}" \
            -i "{input.bam}" &>> "{log}"
        """


rule qc_calculate_three_prime_bias_binned_short_reads:
    input:
        gencode_bed="results/qc/convert_gtf_to_bed_binned/{nam}.bed",
        bam="results/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}.bam",
    output:
        "results/qc/qc_calculate_three_prime_bias_binned_short_reads/{read_number}/{platform}/{sample}/{nam}/three_prime_bias.geneBodyCoverage.txt",
    log:
        "logs/qc/qc_calculate_three_prime_bias_binned_short_reads/{read_number}/{platform}/{sample}/{nam}.log",
    conda:
        "../envs/standalone/rseqc.yaml"
    threads: 1
    params:
        output_path="results/qc/qc_calculate_three_prime_bias_binned_short_reads/{read_number}/{platform}/{sample}/{nam}/three_prime_bias",
    shell:
        """
        conda list > "{log}";
        geneBody_coverage.py \
            -r "{input.gencode_bed}" \
            -o "{params.output_path}" \
            -i "{input.bam}" &>> "{log}"
        """


rule qc_convert_gtf_to_bed:
    input:
        "results/setup/download_transcriptome/gencode.v45.primary_assembly.annotation.gtf",
    output:
        exons="results/qc/qc_convert_gtf_to_bed/exons.bed",
        genes="results/qc/qc_convert_gtf_to_bed/genes.bed",
    log:
        "results/qc/qc_convert_gtf_to_bed/out.log",
    conda:
        "../envs/standalone/read_classes.yaml"
    shell:
        """
        conda list > "{log}"
        awk '$3 == "exon" {{print $1 "\t" $4-1 "\t" $5}}' {input} | sort -k1,1 -k2,2n > {output.exons}
        awk '$3 == "gene" {{print $1 "\t" $4-1 "\t" $5}}' {input} | sort -k1,1 -k2,2n > {output.genes}
        """


rule qc_calculate_read_classes:
    input:
        bam="results/qc/qc_align_genome/{tech}/{read_number}/{platform}/{sample}.bam",
        exons="results/qc/qc_convert_gtf_to_bed/exons.bed",
        genes="results/qc/qc_convert_gtf_to_bed/genes.bed",
    output:
        "results/qc/qc_calculate_read_classes/{tech}/{read_number}/{platform}/{sample}.tsv",
    log:
        "logs/qc/qc_calculate_read_classes/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/read_classes.yaml"
    threads: 24
    shell:
        """
        conda list > "{log}";
        
        samtools view -u -F 2304 {input.bam} | \
            bedtools tag -i stdin -files {input.exons} {input.genes} -labels EX GE -tag GC -f 0.05 | \
            samtools view - | \
            awk 'BEGIN {{OFS="\t"}} {{
                if ($0 ~ /GC:Z:.*EX/) {{ print $1, "Exonic" }} 
                else if ($0 ~ /GC:Z:.*GE/) {{ print $1, "Intronic" }} 
                else {{ print $1, "Intergenic" }}
            }}' > {output}
        """


rule qc_calculate_read_classes_short_reads_:
    input:
        bam="results/qc/qc_align_genome_short_read/{read_number}/{platform}/{sample}.bam",
        exons="results/qc/qc_convert_gtf_to_bed/exons.bed",
        genes="results/qc/qc_convert_gtf_to_bed/genes.bed",
    output:
        tsv="results/qc/qc_calculate_read_classes_short_reads/{read_number}/{platform}/{sample}.tsv",
    log:
        "logs/qc/qc_calculate_read_classes_short_reads/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/standalone/read_classes.yaml"
    threads: 24
    resources:
        mem_gb=8,
    shell:
        """
        conda list > "{log}" 2>&1

        
        if [[ "{wildcards.platform}" == "bulk" ]]; then
            echo "Detected platform 'bulk'. Running Paired-End logic..." >> "{log}"
            
            samtools view -u -F 2304 {input.bam} | \
            bedtools tag -i stdin -files {input.exons} {input.genes} -labels EX GE -tag GC -f 0.5 | \
            samtools view - | \
            awk 'BEGIN {{OFS="\t"}} {{
                if ($0 ~ /GC:Z:.*EX/) {{ cat="Exonic"; score=3 }}
                else if ($0 ~ /GC:Z:.*GE/) {{ cat="Intronic"; score=2 }}
                else {{ cat="Intergenic"; score=1 }}
                print $1, score, cat
            }}' | \
            sort -S 4G -k1,1 -k2,2nr | \
            awk '!seen[$1]++ {{print $1, $3}}' > {output.tsv} 2>> "{log}"

        else
            echo "Detected platform '{wildcards.platform}'..." >> "{log}"

            samtools view -u -F 2304 {input.bam} | \
            bedtools tag -i stdin -files {input.exons} {input.genes} -labels EX GE -tag GC -f 0.5 | \
            samtools view - | \
            awk 'BEGIN {{OFS="\t"}} {{
                if ($0 ~ /GC:Z:.*EX/) {{ print $1, "Exonic" }}
                else if ($0 ~ /GC:Z:.*GE/) {{ print $1, "Intronic" }}
                else {{ print $1, "Intergenic" }}
            }}' > {output.tsv} 2>> "{log}"
        fi
        """


rule qc_calculate_internal_priming:
    input:
        bam="results/qc/qc_align_genome/{tech}/{read_number}/{platform}/{sample}.bam",
        genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        summary="results/qc/qc_calculate_internal_priming/{tech}/{read_number}/{platform}/{sample}_summary.txt",
        gene_count="results/qc/qc_calculate_internal_priming/{tech}/{read_number}/{platform}/{sample}_count.txt",
        bam="results/qc/qc_calculate_internal_priming/{tech}/{read_number}/{platform}/{sample}.bam",
    log:
        "results/qc/qc_calculate_internal_priming/{tech}/{read_number}/{platform}/{sample}.log",
    conda:
        "../envs/py/prime_spotter.yaml"
    threads: 6
    shell:
        """
        conda list > "{log}" 2>&1
        python workflow/scripts/py/PrimeSpotter/PrimeSpotter/PrimeSpotter.py \
            --bam_file {input.bam} \
            --genome-ref {input.genome} \
            --gtf_file {input.transcriptome} \
            --output-summary {output.summary} \
            --output-gene-count {output.gene_count} \
            --processes {threads} | samtools view -@ {threads} -S -b > {output.bam}
        """
