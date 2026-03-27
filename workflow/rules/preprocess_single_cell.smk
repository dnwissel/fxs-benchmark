wildcard_constraints:
    subsample_number=r"[^/]+(?<!full)",


rule demultiplex_sc_ont:
    input:
        "results/downsample_sc/ont/{subsample_number}_{number_to_sample}/{sample}.fastq.gz",
    output:
        "results/preprocess/demultiplex_sc_ont/{subsample_number}_{number_to_sample}/{sample}matched_reads.fastq.gz",
    log:
        "logs/preprocess/demultiplex_sc_ont/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/py/blaze.yaml"
    threads: config["demultiplex_threads"]
    params:
        expect_cells=lambda w: config["demultiplex_expect_cells"][f"{w.sample}"],
        kit_version=config["blaze_sc_kit_version"],
        output_dir="results/preprocess/demultiplex_sc_ont/{subsample_number}_{number_to_sample}/{sample}",
    shell:
        """
        conda list &> {log};
        blaze --expect-cells={params.expect_cells} \
              --output-prefix {params.output_dir} \
              --threads={threads} \
              --overwrite \
              --kit-version {params.kit_version} \
              {input} &>> {log}
        """


rule align_ont_sc_transcriptome:
    input:
        reads="results/preprocess/demultiplex_sc_ont/{subsample_number}_{number_to_sample}/{sample}matched_reads.fastq.gz",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
    output:
        "results/preprocess/align_ont_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/preprocess/align_ont_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: config["align_map_bam_threads"]
    params:
        map_threads=config["align_map_bam_threads"],
        sort_threads=config["align_sort_bam_threads"],
        sort_memory_gb=config["align_sort_bam_memory_gb"],
    shell:
        """
        minimap2 -y --eqx -N 100 -ax map-ont \
            -t {params.map_threads} \
            {input.gencode_transcriptome} {input.reads} \
            | samtools sort \
            -@ {params.sort_threads} \
            -m{params.sort_memory_gb}g \
            -o {output} --write-index - > {log} \
            2> {log}
        """


rule filter_primary_alignments:
    input:
        bam="results/preprocess/align_ont_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    output:
        bam="results/preprocess/filter_primary_alignments/{subsample_number}_{number_to_sample}/{sample}.primary.bam",
        bai="results/preprocess/filter_primary_alignments/{subsample_number}_{number_to_sample}/{sample}.primary.bam.bai",
    log:
        "logs/preprocess/filter_primary_alignments/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["align_map_bam_threads"]
    params:
        map_threads=config["align_map_bam_threads"],
        sort_threads=config["align_sort_bam_threads"],
        sort_memory_gb=config["align_sort_bam_memory_gb"],
    shell:
        """
        conda list &> {log};
        samtools view \
            -@ {threads} \
            -F 2304 \
            -b \
            -o {output.bam} \
            {input.bam}

        samtools index \
            -@ {threads} \
            {output.bam} > {log} 2>&1
        """


rule prepare_gene_transcript_map_umi_tools:
    input:
        "results/setup/standardize_gtf_files/gencode_map_headered.txt",
    output:
        "results/preprocess/prepare_gene_transcript_map_umi_tools/umi_tools_gene_transcript_map.tsv",
    log:
        "logs/preprocess/prepare_gene_transcript_map_umi_tools/umi_tools_gene_transcript_map.log",
    conda:
        "../envs/standalone/seqtk.yaml"
    shell:
        """
        awk 'BEGIN{{OFS="\t"}} NR > 1 {{print $2, $1}}' {input} \
            > {output} 2>> {log}
        """


rule deduplicate_run_umi_tools_transcriptome:
    input:
        reads="results/preprocess/filter_primary_alignments/{subsample_number}_{number_to_sample}/{sample}.primary.bam",
        gene_transcript_map="results/preprocess/prepare_gene_transcript_map_umi_tools/umi_tools_gene_transcript_map.tsv",
    output:
        "results/preprocess/deduplicate_run_umi_tools_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/preprocess/deduplicate_run_umi_tools_transcriptome/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/umi_tools.yaml"
    threads: 8
    params:
        seed=config["seed"],
    shell:
        """
        conda list &> {log};
        export PYTHONHASHSEED={params.seed};
        umi_tools dedup \
            --per-gene \
            --per-contig \
            --gene-transcript-map={input.gene_transcript_map} \
            --stdin={input.reads} \
            --log={log} \
            --per-cell \
            --extract-umi-method=tag \
            --umi-tag=UB \
            --random-seed {params.seed} \
            --cell-tag CB > {output} 2>> {log}
        """


rule restore_secondary_alignments:
    input:
        dedup_bam="results/preprocess/deduplicate_run_umi_tools_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
        original_bam="results/preprocess/align_ont_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    output:
        bam="results/preprocess/restore_secondary_alignments/{subsample_number}_{number_to_sample}/{sample}.dedup_with_secondary.bam",
        bai="results/preprocess/restore_secondary_alignments/{subsample_number}_{number_to_sample}/{sample}.dedup_with_secondary.bam.bai",
    log:
        "logs/preprocess/restore_secondary_alignments/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: config["stall_io_threads"]
    params:
        read_ids=lambda wildcards: f"results/preprocess/restore_secondary_alignments/{wildcards.subsample_number}_{wildcards.number_to_sample}/{wildcards.sample}.read_ids.txt",
    shell:
        """
        samtools view -@ {threads} {input.dedup_bam} | cut -f 1 | sort -u -T $(dirname {output.bam}) > {params.read_ids} 2> {log}

        samtools view \
            -@ {threads} \
            -N {params.read_ids} \
            -o {output.bam} \
            -b \
            {input.original_bam} 2>> {log}

        samtools index -@ {threads} {output.bam} 2>> {log}

        rm {params.read_ids}
        """


rule preprocess_pb:
    input:
        reads="results/downsample_sc/pb/{subsample_number}_{number_to_sample}/{sample}.bam",
        primers="config/pb_3p_primers.fasta",
        barcodes="config/3M-february-2018-REVERSE-COMPLEMENTED.txt",
    output:
        "results/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.deduplicated.bam",
    log:
        "logs/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/isoseq.yaml"
    threads: config["stall_io_threads"]
    params:
        outfile_primers="results/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.bam",
        primers="results/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.5p--3p.bam",
        tag="results/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.fltn.bam",
        refine="results/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.fltnc.bam",
        correct="results/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.fltnc.corrected.bam",
        sort_threads=config["align_sort_bam_threads"],
        sort_memory_gb=config["align_sort_bam_memory_gb"],
        sort="results/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.fltnc.corrected.sorted.bam",
    shell:
        """
        conda list &> {log};
        lima {input.reads} {input.primers} {params.outfile_primers} --no-reports \
            --num-threads {threads} --isoseq &>> {log};
        isoseq tag {params.primers} {params.tag} \
            --design T-12U-16B --num-threads {threads} &>> {log};
        rm {params.primers} &>> {log};
        isoseq refine {params.tag} {input.primers} {params.refine} \
            --require-polya --num-threads {threads};
        rm {params.tag} &>> {log};
        isoseq correct --barcodes {input.barcodes} \
            --num-threads {threads} {params.refine} {params.correct} \
            &>> {log};
        rm {params.refine} &>> {log};
        samtools sort -m{params.sort_memory_gb}g \
            -@ {params.sort_threads}  -t CB {params.correct} \
             -o {params.sort} &>> {log};
        rm {params.correct} &>> {log};
        isoseq groupdedup {params.sort} {output} &>> {log};
        rm {params.sort} &>> {log}
        """


rule align_pb_sc_transcriptome:
    input:
        reads="results/downsample_sc/preprocess_pb/{subsample_number}_{number_to_sample}/{sample}.deduplicated.bam",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
    output:
        "results/preprocess/align_pb_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/preprocess/align_pb_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: 24
    params:
        map_threads=config["align_map_bam_threads"],
        sort_threads=config["align_sort_bam_threads"],
        sort_memory_gb=config["align_sort_bam_memory_gb"],
    shell:
        """
        samtools fastq -T CB -@ {threads} {input.reads} | \
            minimap2 -y --eqx -N 100 -ax map-hifi \
            -t {params.map_threads} \
            {input.gencode_transcriptome} - \
            | samtools sort \
            -@ {params.sort_threads} \
            -m{params.sort_memory_gb}g \
            -o {output} - > {log} \
            2> {log}
        """


rule align_pb_sc_genome:
    input:
        reads="results/preprocess/align_pb_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_bed="results/setup/convert_gtfs_to_beds/gencode.bed",
    output:
        "results/preprocess/align_pb_sc_genome/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/preprocess/align_pb_sc_genome/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: 24
    params:
        map_threads=config["align_map_bam_threads"],
        sort_threads=config["align_sort_bam_threads"],
        sort_memory_gb=config["align_sort_bam_memory_gb"],
    shell:
        """
        samtools fastq -T CB -@ {threads} {input.reads} | \
            minimap2 -y -ax splice:hq -uf \
            --junc-bed {input.gencode_bed} \
            -t {params.map_threads} \
            {input.gencode_genome} - \
            | samtools sort \
            -@ {params.sort_threads} \
            -m{params.sort_memory_gb}g \
            --write-index \
            -o {output} - > {log} \
            2> {log}
        """


rule align_ont_sc_genome:
    input:
        reads="results/preprocess/restore_secondary_alignments/{subsample_number}_{number_to_sample}/{sample}.dedup_with_secondary.bam",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_bed="results/setup/convert_gtfs_to_beds/gencode.bed",
    output:
        "results/preprocess/align_ont_sc_genome/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/preprocess/align_ont_sc_genome/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: 24
    params:
        map_threads=config["align_map_bam_threads"],
        sort_threads=config["align_sort_bam_threads"],
        sort_memory_gb=config["align_sort_bam_memory_gb"],
    shell:
        """
        samtools fastq -T CB -@ {threads} {input.reads} | \
            minimap2 -y -ax splice -uf \
            --junc-bed {input.gencode_bed} \
            -t {params.map_threads} \
            {input.gencode_genome} - \
            | samtools sort \
            -@ {params.sort_threads} \
            -m{params.sort_memory_gb}g \
            --write-index \
            -o {output} - > {log} \
            2> {log}
        """


rule prepare_kallisto_sc_pb_fastq:
    input:
        "results/preprocess/align_pb_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    output:
        "results/preprocess/prepare_kallisto_sc_pb_fastq/{subsample_number}_{number_to_sample}/{sample}.fastq.gz",
    log:
        "logs/preprocess/prepare_kallisto_sc_pb_fastq/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/pigz.yaml"
    threads: 12
    shell:
        """
        samtools fastq -@ 4 -T CB {input} | \
        awk '{{
            if (NR%4==1) {{
                n = split($0, a, " ");
                header = a[1];
                split(a[2], b, ":");
                barcode = b[3];
                print header;
            }}
            if (NR%4==2) {{
                print barcode $0;
            }}
            if (NR%4==3) {{
                print $0;
            }}
            if (NR%4==0) {{
                qual_placeholder = sprintf("%*s", length(barcode), "");
                gsub(/ /, "I", qual_placeholder);
                print qual_placeholder $0;
            }}
        }}' | pigz -p 4 > {output} 2>> {log}
        """


rule prepare_kallisto_sc_ont_fastq:
    input:
        "results/preprocess/restore_secondary_alignments/{subsample_number}_{number_to_sample}/{sample}.dedup_with_secondary.bam",
    output:
        "results/preprocess/prepare_kallisto_sc_ont_fastq/{subsample_number}_{number_to_sample}/{sample}.fastq.gz",
    log:
        "logs/preprocess/prepare_kallisto_sc_ont_fastq/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/pigz.yaml"
    threads: 12
    shell:
        """
        samtools fastq -@ 4 -T CB {input} | \
        awk '{{
            if (NR%4==1) {{
                n = split($0, a, " ");
                header = a[1];
                split(a[2], b, ":");
                barcode = b[3];
                print header;
            }}
            if (NR%4==2) {{
                print barcode $0;
            }}
            if (NR%4==3) {{
                print $0;
            }}
            if (NR%4==0) {{
                qual_placeholder = sprintf("%*s", length(barcode), "");
                gsub(/ /, "I", qual_placeholder);
                print qual_placeholder $0;
            }}
        }}' | pigz -p 4 > {output} 2>> {log}
        """


rule sort_ont_sc_transcriptome:
    input:
        "results/preprocess/restore_secondary_alignments/{subsample_number}_{number_to_sample}/{sample}.dedup_with_secondary.bam",
    output:
        "results/preprocess/sort_ont_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/preprocess/sort_ont_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: 24
    params:
        map_threads=config["align_map_bam_threads"],
        sort_threads=config["align_sort_bam_threads"],
        sort_memory_gb=config["align_sort_bam_memory_gb"],
    shell:
        """
        conda list &> {log};
        samtools sort -t CB -n {input} -o {output} \
            -m{params.sort_memory_gb}g -@{params.sort_threads} &>> {log}
        """


rule sort_pb_sc_transcriptome:
    input:
        "results/preprocess/align_pb_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    output:
        "results/preprocess/sort_pb_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/preprocess/sort_pb_sc_transcriptome/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: 24
    params:
        map_threads=config["align_map_bam_threads"],
        sort_threads=config["align_sort_bam_threads"],
        sort_memory_gb=config["align_sort_bam_memory_gb"],
    shell:
        """
        conda list &> {log};
        samtools sort -t CB -n {input} -o {output} \
            -m{params.sort_memory_gb}g -@{params.sort_threads} &>> {log}
        """


rule adjust_lr_sc_genome:
    input:
        "results/preprocess/align_{tech}_sc_genome/{subsample_number}_{number_to_sample}/{sample}.bam",
    output:
        "results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam",
    log:
        "logs/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/jvar.yaml"
    threads: 24
    shell:
        """
        conda list &> {log};
        samjdk -e \
            'final Object s=record.getAttribute(\\"CB\\"); if(s!=null) {{record.setAttribute(\\"CB\\",null);record.setAttribute(\\"BC\\",s);}} return record;' \
            {input} -o {output} &>> {log}
        """


rule index_adjust_lr_sc_genome:
    input:
        "results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam",
    output:
        "results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam.bai",
    log:
        "logs/preprocess/index_adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/samtools.yaml"
    threads: 24
    shell:
        """
        conda list &> {log};
        samtools index {input} &>> {log}
        """
