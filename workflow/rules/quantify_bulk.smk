configfile: "config/config.yaml"


rule quantify_bulk_joint_run_miniquant_lr_ont:
    input:
        ont_reads=f"{config["bulk_ont_directory"]}/{{sample}}_R1.fastq.gz",
        illumina_reads_first=f"{config["bulk_ill_directory"]}/{{sample}}_R1.fastq.gz",
        illumina_reads_second=f"{config["bulk_ill_directory"]}/{{sample}}_R2.fastq.gz",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
    output:
        "results/quantify_bulk_joint/ont/run_miniquant_lr/{sample}/abundance.tsv",
        directory("results/quantify_bulk_joint/ont/run_miniquant_lr/{sample}"),
    log:
        "logs/quantify_bulk_joint/ont/run_miniquant_lr/{sample}.log",
    container:
        "docker://tidesun/miniquant@sha256:69c2aadb580263beb0d236649d4c3a65100b617b79eb7cd85e6c4ed2e3e8931e"
    threads: 12
    params:
        outdir="results/quantify_bulk_joint/ont/run_miniquant_lr/{sample}",
    shell:
        """
        /app/miniQuant_linux/miniQuant -h &> {log} || true;
        /app/miniQuant_linux/miniQuant quant -r {input.gencode_transcriptome} \
            -l {input.ont_reads} -t {threads} \
            -1 {input.illumina_reads_first} \
            -2 {input.illumina_reads_second} \
            -o {params.outdir} \
            --long_reads_library_prep cDNA-ONT \
            --short_reads_strandness rf-stranded \
            &>> {log}
        """


rule quantify_bulk_joint_run_miniquant_lr_pb:
    input:
        pb_reads="results/prepare_bulk/pb_convert_to_bam/{sample}.fastq.gz",
        illumina_reads_first=f"{config["bulk_ill_directory"]}/{{sample}}_R1.fastq.gz",
        illumina_reads_second=f"{config["bulk_ill_directory"]}/{{sample}}_R2.fastq.gz",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
    output:
        "results/quantify_bulk_joint/pb/run_miniquant_lr/{sample}/abundance.tsv",
        directory("results/quantify_bulk_joint/pb/run_miniquant_lr/{sample}"),
    log:
        "logs/quantify_bulk_joint/pb/run_miniquant_lr/{sample}.log",
    container:
        "docker://tidesun/miniquant@sha256:69c2aadb580263beb0d236649d4c3a65100b617b79eb7cd85e6c4ed2e3e8931e"
    threads: 12
    params:
        outdir="results/quantify_bulk_joint/pb/run_miniquant_lr/{sample}",
    shell:
        """
        /app/miniQuant_linux/miniQuant -h &> {log} || true;
        /app/miniQuant_linux/miniQuant quant -r {input.gencode_transcriptome} \
            -l {input.pb_reads} -t {threads} \
            -1 {input.illumina_reads_first} \
            -2 {input.illumina_reads_second} \
            -o {params.outdir} \
            --long_reads_library_prep cDNA-PacBio \
            --short_reads_strandness rf-stranded \
            &>> {log}
        """


rule quantify_bulk_downsampled_run_salmon_illumina:
    input:
        gencode_index="results/setup/create_salmon_index/gencode",
        sirv_index="results/setup/create_salmon_index/sirv",
        first_reads="results/downsample_bulk/short_read_fastq/{subsample_number}_{read_number}_{data_type}/{sample}-r1.fastq.gz",
        second_reads="results/downsample_bulk/short_read_fastq/{subsample_number}_{read_number}_{data_type}/{sample}-r2.fastq.gz",
    output:
        "results/quantify_bulk_downsampled/run_salmon_illumina/{subsample_number}_{read_number}/{data_type}/{sample}/illumina/quant.sf",
        directory(
            "results/quantify_bulk_downsampled/run_salmon_illumina/{subsample_number}_{read_number}/{data_type}/{sample}/illumina"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_salmon_illumina/{subsample_number}_{read_number}_{data_type}/{sample}.log",
    conda:
        "../envs/standalone/salmon.yaml"
    threads: config["quantify_isoforms_threads"]
    params:
        output_name="results/quantify_bulk_downsampled/run_salmon_illumina/{subsample_number}_{read_number}/{data_type}/{sample}/illumina",
        fragment_length_mean=config["salmon_fragment_length_mean"],
        fragment_length_sd=config["salmon_fragment_length_sd"],
        num_bootstraps=config["num_bootstraps"],
    shell:
        """
        conda list > {log};
        if [[ "{wildcards.data_type}" == "gencode" ]]; then
            salmon quant -i {input.gencode_index} -l A \
                -1 {input.first_reads} -2 {input.second_reads} \
                --validateMappings -o {params.output_name} -p {threads} \
                --seqBias --gcBias --fldMean {params.fragment_length_mean} \
                --numBootstraps {params.num_bootstraps} \
                --fldSD {params.fragment_length_sd} &>> {log}
        else
            salmon quant -i {input.sirv_index} -l A \
                -1 {input.first_reads} -2 {input.second_reads} \
                --validateMappings -o {params.output_name} -p {threads} \
                --seqBias --gcBias --fldMean {params.fragment_length_mean} \
                --numBootstraps {params.num_bootstraps} \
                --fldSD {params.fragment_length_sd} &>> {log}
        fi
        """


rule quantify_bulk_downsampled_dump_salmon_corrected:
    input:
        "results/quantify_bulk_downsampled/run_salmon_illumina/{subsample_number}_{read_number}/{data_type}/{sample}/illumina",
    output:
        "results/quantify_bulk_downsampled/run_salmon_illumina_corrected/{subsample_number}_{read_number}/{data_type}/{sample}/illumina/quant.sf",
        directory(
            "results/quantify_bulk_downsampled/run_salmon_illumina_corrected/{subsample_number}_{read_number}/{data_type}/{sample}/illumina"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_salmon_illumina_corrected/{subsample_number}_{read_number}/{data_type}/{sample}/illumina/{sample}.log",
    conda:
        "../envs/r/salmon.yaml"
    threads: 1
    script:
        "../scripts/r/dump_salmon_illumina.R"


rule quantify_bulk_downsampled_align_genome_bulk_lr:
    input:
        reads="results/downsample_bulk/long_read_fastq/{type}/{subsample_number}_{read_number}_{data_type}/{sample}-r1.fastq.gz",
        sirv_genome="results/setup/adjust_sirv_names/sirv_set_four.fa",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_bed="results/setup/convert_gtfs_to_beds/gencode.bed",
        sirv_bed="results/setup/convert_gtfs_to_beds/sirv.bed",
    output:
        "results/quantify_downsampled/align_genome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.aligned.sorted.bam",
    log:
        "logs/quantify_downsampled/align_genome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: config["align_threads"]
    params:
        align_sort_bam_memory_gb=config["align_sort_bam_memory_gb"],
        align_sort_bam_threads=config["align_sort_bam_threads"],
    shell:
        """
        conda list &> {log};
        if [[ "{wildcards.data_type}" != "gencode" && "{wildcards.type}" == "pb" ]]; then
            minimap2 -ax splice:hq --junc-bed {input.sirv_bed} \
                --splice-flank=no \
                -t {threads} -uf {input.sirv_genome} \
                {input.reads} 2>> {log} | \
                samtools sort -@ {params.align_sort_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g -o {output} \
                - &>> {log};
            sleep 61s &>> {log};
            samtools index {output} &>> {log}
        elif [[ "{wildcards.data_type}" == "gencode" && "{wildcards.type}" == "pb" ]]; then
            minimap2 -ax splice:hq --junc-bed {input.gencode_bed} \
                -t {threads} -uf {input.gencode_genome} \
                {input.reads} 2>> {log} | \
                samtools sort -@ {params.align_sort_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g -o {output} \
                - &>> {log};
            sleep 61s &>> {log};
            samtools index {output} &>> {log}
        elif [[ "{wildcards.data_type}" == "gencode" && "{wildcards.type}" == "ont" ]]; then
            minimap2 -ax splice --junc-bed {input.gencode_bed} \
                -t {threads} -uf {input.gencode_genome} \
                {input.reads} 2>> {log} | \
                samtools sort -@ {params.align_sort_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g -o {output} \
                - &>> {log};
            sleep 61s &>> {log};
            samtools index {output} &>> {log}
        else
            minimap2 -ax splice --junc-bed {input.sirv_bed} \
                --splice-flank=no \
                -t {threads} -uf {input.sirv_genome} \
                {input.reads} 2>> {log} | \
                samtools sort -@ {params.align_sort_bam_threads} \
                -m{params.align_sort_bam_memory_gb}g -o {output} \
                - &>> {log};
            sleep 61s &>> {log};
            samtools index {output} &>> {log}
        fi
        """


rule quantify_bulk_downsampled_align_transcriptome_bulk_lr:
    input:
        reads="results/downsample_bulk/long_read_fastq/{type}/{subsample_number}_{read_number}_{data_type}/{sample}-r1.fastq.gz",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
        sirv_transcriptome="results/setup/extract_transcriptomes/sirv_transcriptome.fa",
    output:
        "results/quantify_downsampled/align_transcriptome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.aligned.bam",
    log:
        "logs/quantify_downsampled/align_transcriptome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    threads: config["align_threads"]
    params:
        align_sort_bam_threads=config["align_sort_bam_threads"],
        n_secondary_alignments=100,
    shell:
        """
        conda list &> {log};
        if [[ "{wildcards.data_type}" != "gencode" && "{wildcards.type}" == "pb" ]]; then
            minimap2 --eqx -N {params.n_secondary_alignments} -ax map-hifi \
                -t {threads} {input.sirv_transcriptome} \
                {input.reads} 2>> {log} | \
                samtools view -@ {params.align_sort_bam_threads} \
                -bo {output} \
                - &>> {log};
        elif [[ "{wildcards.data_type}" == "gencode" && "{wildcards.type}" == "pb" ]]; then
            minimap2 --eqx -N {params.n_secondary_alignments} -ax map-hifi \
                -t {threads} {input.gencode_transcriptome} \
                {input.reads} 2>> {log} | \
                samtools view -@ {params.align_sort_bam_threads} \
                -bo {output} \
                - &>> {log};
        elif [[ "{wildcards.data_type}" == "gencode" && "{wildcards.type}" == "ont" ]]; then
            minimap2 --eqx -N {params.n_secondary_alignments} -ax lr:hq \
                -t {threads} {input.gencode_transcriptome} \
                {input.reads} 2>> {log} | \
                samtools view -@ {params.align_sort_bam_threads} \
                -bo {output} \
                - &>> {log};
        else
            minimap2 --eqx -N {params.n_secondary_alignments} -ax lr:hq \
                -t {threads} {input.sirv_transcriptome} \
                {input.reads} 2>> {log} | \
                samtools view -@ {params.align_sort_bam_threads} \
                -bo {output} \
                - &>> {log};
        fi
        """


rule quantify_bulk_downsampled_run_bambu_lr:
    input:
        "results/setup/install_bambu/done.txt",
        reads="results/quantify_downsampled/align_genome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.aligned.sorted.bam",
        sirv_genome="results/setup/adjust_sirv_names/sirv_set_four.fa",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/setup/adjust_sirv_names/sirv_set_four.gtf",
    output:
        "results/quantify_bulk_downsampled/run_bambu_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/counts_transcript.txt",
        directory(
            "results/quantify_bulk_downsampled/run_bambu_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_bambu_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    benchmark:
        "benchmarks/quantify_bulk_downsampled/run_bambu_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log"
    conda:
        "../envs/r/bambu.yaml"
    threads: config["quantify_isoforms_threads"]
    params:
        discovery=False,
        NDR=None,
        quantification=True,
        outdir="results/quantify_bulk_downsampled/run_bambu_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}",
    script:
        "../scripts/r/run_bambu.R"


rule quantify_bulk_downsampled_run_isoquant_lr:
    input:
        "results/setup/install_bambu/done.txt",
        reads="results/quantify_downsampled/align_genome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.aligned.sorted.bam",
        sirv_genome="results/setup/adjust_sirv_names/sirv_set_four.fa",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/setup/adjust_sirv_names/sirv_set_four.gtf",
    output:
        "results/quantify_bulk_downsampled/run_isoquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/OUT/OUT.transcript_counts.tsv",
        directory(
            "results/quantify_bulk_downsampled/run_isoquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_isoquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    benchmark:
        "benchmarks/quantify_bulk_downsampled/run_isoquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log"
    conda:
        "../envs/py/isoquant.yaml"
    threads: config["quantify_isoforms_threads"]
    params:
        discovery=False,
        NDR=None,
        quantification=True,
        outdir="results/quantify_bulk_downsampled/run_isoquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}",
    shell:
        """
        conda list &> {log};
        if [[ "{wildcards.data_type}" != "gencode" ]]; then
            if [[ "{wildcards.type}" == "ont" ]]; then
                isoquant.py --reference {input.sirv_genome} \
                    --genedb {input.sirv_transcriptome} \
                    --bam {input.reads} \
                    --threads {threads} \
                    --data_type ont \
                    --polya_requirement never \
                    --no_model_construction \
                    -o {params.outdir} &>> {log}
            else
                isoquant.py --reference {input.sirv_genome} \
                    --genedb {input.sirv_transcriptome} \
                    --bam {input.reads} \
                    --threads {threads} \
                    --data_type pacbio \
                    --polya_requirement never \
                    --no_model_construction \
                    -o {params.outdir} &>> {log}
            fi
        else
          if [[ "{wildcards.type}" == "ont" ]]; then
                isoquant.py --reference {input.gencode_genome} \
                    --genedb {input.gencode_transcriptome} \
                    --complete_genedb \
                    --bam {input.reads} \
                    --threads {threads} \
                    --data_type ont \
                    --polya_requirement never \
                    --no_model_construction \
                    -o {params.outdir} &>> {log}
            else
                isoquant.py --reference {input.gencode_genome} \
                    --genedb {input.gencode_transcriptome} \
                    --complete_genedb \
                    --bam {input.reads} \
                    --threads {threads} \
                    --data_type pacbio \
                    --polya_requirement never \
                    --no_model_construction \
                    -o {params.outdir} &>> {log}
            fi
        fi
        """


rule quantify_bulk_downsampled_run_isosceles_lr:
    input:
        "results/setup/install_isosceles/done.txt",
        reads="results/quantify_downsampled/align_genome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.aligned.sorted.bam",
        sirv_genome="results/setup/adjust_sirv_names/sirv_set_four.fa",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
        sirv_transcriptome="results/setup/adjust_sirv_names/sirv_set_four.gtf",
    output:
        "results/quantify_bulk_downsampled/run_isosceles_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/counts.tsv",
        directory(
            "results/quantify_bulk_downsampled/run_isosceles_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_isosceles_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    benchmark:
        "benchmarks/quantify_bulk_downsampled/run_isosceles_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log"
    conda:
        "../envs/r/isosceles.yaml"
    threads: config["quantify_isoforms_threads"]
    params:
        seed=config["seed"],
    script:
        "../scripts/r/run_isosceles.R"


rule quantify_bulk_downsampled_run_kallisto_long_lr:
    input:
        reads="results/downsample_bulk/long_read_fastq/{type}/{subsample_number}_{read_number}_{data_type}/{sample}-r1.fastq.gz",
        sirv_idx=f"results/setup/create_lr_kallisto_index/sirv_k-{config['lr_kallisto_index_k']}.idx",
        gencode_idx=f"results/setup/create_lr_kallisto_index/gencode_k-{config['lr_kallisto_index_k']}.idx",
        gencode_transcriptome_gmap_headered="results/setup/standardize_gtf_files/gencode_map_headered.txt",
        sirv_transcriptome_gmap_headered="results/setup/standardize_gtf_files/sirv_set_four_map_headered.txt",
    output:
        "results/quantify_bulk_downsampled/run_kallisto_long_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/abundance_1.tsv",
        directory(
            "results/quantify_bulk_downsampled/run_kallisto_long_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_kallisto_long_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    benchmark:
        "benchmarks/quantify_bulk_downsampled/run_kallisto_long_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log"
    conda:
        "../envs/standalone/kallisto.yaml"
    threads: config["quantify_isoforms_threads"]
    params:
        bus_threshold=config["kallisto_bus_threshold"],
        seed=config["seed"],
        outdir="results/quantify_bulk_downsampled/run_kallisto_long_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}",
    shell:
        """
        conda list &> {log};
        if [[ "{wildcards.data_type}" != "gencode" ]]; then
            if [[ "{wildcards.type}" == "ont" ]]; then
                kallisto bus -t {threads} --long --threshold \
                    {params.bus_threshold} -x bulk -i {input.sirv_idx} \
                    -o {params.outdir} {input.reads} &>> {log};
                bustools sort -t {threads} {params.outdir}/output.bus \
                    -o {params.outdir}/sorted.bus &>> {log};
                bustools count {params.outdir}/sorted.bus \
                    -t {params.outdir}/transcripts.txt \
                    -e {params.outdir}/matrix.ec \
                    -g {input.sirv_transcriptome_gmap_headered} \
                    -o {params.outdir}/count --cm -m \
                    &>> {log};
                kallisto quant-tcc -t {threads} \
                    --long -P ONT -f {params.outdir}/flens.txt \
                    {params.outdir}/count.mtx -i {input.sirv_idx} \
                    -e {params.outdir}/count.ec.txt \
                    -o {params.outdir} \
                    --seed {params.seed} \
                    --matrix-to-files \
                    &>> {log}
            else
                kallisto bus -t {threads} --long --threshold \
                            {params.bus_threshold} -x bulk -i {input.sirv_idx} \
                            -o {params.outdir} {input.reads} &>> {log};
                        bustools sort -t {threads} {params.outdir}/output.bus \
                            -o {params.outdir}/sorted.bus &>> {log};
                        bustools count {params.outdir}/sorted.bus \
                            -t {params.outdir}/transcripts.txt \
                            -e {params.outdir}/matrix.ec \
                            -g {input.sirv_transcriptome_gmap_headered} \
                            -o {params.outdir}/count --cm -m \
                            &>> {log};
                        kallisto quant-tcc -t {threads} \
                            --long -P PacBio -f {params.outdir}/flens.txt \
                            {params.outdir}/count.mtx -i {input.sirv_idx} \
                            -e {params.outdir}/count.ec.txt \
                            -o {params.outdir} \
                            --seed {params.seed} \
                            --matrix-to-files \
                            &>> {log}
            fi
        else
            if [[ "{wildcards.type}" == "ont" ]]; then
                kallisto bus -t {threads} --long --threshold \
                    {params.bus_threshold} -x bulk -i {input.gencode_idx} \
                    -o {params.outdir} {input.reads} &>> {log};
                bustools sort -t {threads} {params.outdir}/output.bus \
                    -o {params.outdir}/sorted.bus &>> {log};
                bustools count {params.outdir}/sorted.bus \
                    -t {params.outdir}/transcripts.txt \
                    -e {params.outdir}/matrix.ec \
                    -g {input.gencode_transcriptome_gmap_headered} \
                    -o {params.outdir}/count --cm -m \
                    &>> {log};
                kallisto quant-tcc -t {threads} \
                    --long -P ONT -f {params.outdir}/flens.txt \
                    {params.outdir}/count.mtx -i {input.gencode_idx} \
                    -e {params.outdir}/count.ec.txt \
                    -o {params.outdir} \
                    --seed {params.seed} \
                    --matrix-to-files \
                    &>> {log}
            else
                kallisto bus -t {threads} --long --threshold \
                    {params.bus_threshold} -x bulk -i {input.gencode_idx} \
                    -o {params.outdir} {input.reads} &>> {log};
                bustools sort -t {threads} {params.outdir}/output.bus \
                    -o {params.outdir}/sorted.bus &>> {log};
                bustools count {params.outdir}/sorted.bus \
                    -t {params.outdir}/transcripts.txt \
                    -e {params.outdir}/matrix.ec \
                    -g {input.gencode_transcriptome_gmap_headered} \
                    -o {params.outdir}/count --cm -m \
                    &>> {log};
                kallisto quant-tcc -t {threads} \
                    --long -P PacBio -f {params.outdir}/flens.txt \
                    {params.outdir}/count.mtx -i {input.gencode_idx} \
                    -e {params.outdir}/count.ec.txt \
                    -o {params.outdir} \
                    --seed {params.seed} \
                    --matrix-to-files \
                    &>> {log}
            fi
        fi
        """


rule quantify_bulk_downsampled_run_miniquant_lr:
    input:
        reads="results/downsample_bulk/long_read_fastq/{type}/{subsample_number}_{read_number}_{data_type}/{sample}-r1.fastq.gz",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
        sirv_transcriptome="results/setup/extract_transcriptomes/sirv_transcriptome.fa",
    output:
        "results/quantify_bulk_downsampled/run_miniquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/abundance.tsv",
        directory(
            "results/quantify_bulk_downsampled/run_miniquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_miniquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    benchmark:
        "benchmarks/quantify_bulk_downsampled/run_miniquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log"
    container:
        "docker://tidesun/miniquant@sha256:69c2aadb580263beb0d236649d4c3a65100b617b79eb7cd85e6c4ed2e3e8931e"
    threads: config["quantify_isoforms_threads"]
    params:
        bus_threshold=config["kallisto_bus_threshold"],
        seed=config["seed"],
        outdir="results/quantify_bulk_downsampled/run_miniquant_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}",
    shell:
        """
        /app/miniQuant_linux/miniQuant -h &> {log} || true;
        if [[ "{wildcards.data_type}" != "gencode" ]]; then
            if [[ "{wildcards.type}" == "ont" ]]; then
                /app/miniQuant_linux/miniQuant quant -r {input.sirv_transcriptome} \
                    -l {input.reads} -t {threads} \
                    -o {params.outdir} \
                    --long_reads_library_prep cDNA-ONT \
                    &>> {log}
            else
                /app/miniQuant_linux/miniQuant quant -r {input.sirv_transcriptome} \
                    -l {input.reads} -t {threads} \
                    -o {params.outdir} \
                    --long_reads_library_prep cDNA-PacBio \
                    &>> {log}
            fi
        else
            if [[ "{wildcards.type}" == "ont" ]]; then
                /app/miniQuant_linux/miniQuant quant -r {input.gencode_transcriptome} \
                    -l {input.reads} -t {threads} \
                    -o {params.outdir} \
                    --long_reads_library_prep cDNA-ONT \
                    &>> {log}
            else
                /app/miniQuant_linux/miniQuant quant -r {input.gencode_transcriptome} \
                    -l {input.reads} -t {threads} \
                    -o {params.outdir} \
                    --long_reads_library_prep cDNA-PacBio \
                    &>> {log}
            fi
        fi
        """


rule quantify_bulk_downsampled_run_oarfish_lr:
    input:
        reads="results/quantify_downsampled/align_transcriptome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.aligned.bam",
    output:
        "results/quantify_bulk_downsampled/run_oarfish_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.quant",
        directory(
            "results/quantify_bulk_downsampled/run_oarfish_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_oarfish_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    benchmark:
        "benchmarks/quantify_bulk_downsampled/run_oarfish_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log"
    conda:
        "../envs/standalone/oarfish.yaml"
    threads: config["quantify_isoforms_threads"]
    params:
        seed=config["seed"],
        bin_width=100,
        filter_group="no-filters",
        outdir="results/quantify_bulk_downsampled/run_oarfish_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}",
    shell:
        """
        conda list &> {log};
        oarfish --threads {threads} \
            --filter-group {params.filter_group} \
            --model-coverage \
            --bin-width {params.bin_width} \
            --alignments {input.reads} \
            --output {params.outdir}/{wildcards.sample} 2> {log}
        """


rule quantify_bulk_downsampled_run_oarfish_lr_bootstrapped:
    input:
        reads="results/quantify_downsampled/align_transcriptome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.aligned.bam",
    output:
        "results/quantify_bulk_downsampled/run_oarfish_lr_bootstrapped/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.quant",
        directory(
            "results/quantify_bulk_downsampled/run_oarfish_lr_bootstrapped/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_oarfish_lr_bootstrapped/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    conda:
        "../envs/standalone/oarfish.yaml"
    threads: 24
    params:
        seed=config["seed"],
        bin_width=100,
        n_bootstraps=30,
        filter_group="no-filters",
        outdir="results/quantify_bulk_downsampled/run_oarfish_lr_bootstrapped/{subsample_number}_{read_number}/{data_type}/{sample}/{type}",
    shell:
        """
        conda list &> {log};
        oarfish --threads {threads} \
            --filter-group {params.filter_group} \
            --model-coverage \
            --num-bootstraps {params.n_bootstraps} \
            --bin-width {params.bin_width} \
            --alignments {input.reads} \
            --output {params.outdir}/{wildcards.sample} 2> {log}
        """


rule quantify_bulk_downsampled_run_salmon_lr:
    input:
        reads="results/quantify_downsampled/align_transcriptome_bulk_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.aligned.bam",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
        sirv_transcriptome="results/setup/extract_transcriptomes/sirv_transcriptome.fa",
    output:
        "results/quantify_bulk_downsampled/run_salmon_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/quant.sf",
        directory(
            "results/quantify_bulk_downsampled/run_salmon_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
        ),
    log:
        "logs/quantify_bulk_downsampled/run_salmon_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log",
    benchmark:
        "benchmarks/quantify_bulk_downsampled/run_salmon_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}/{sample}.log"
    conda:
        "../envs/standalone/salmon.yaml"
    threads: config["quantify_isoforms_threads"]
    params:
        seed=config["seed"],
        bin_width=100,
        filter_group="no-filters",
        outdir="results/quantify_bulk_downsampled/run_salmon_lr/{subsample_number}_{read_number}/{data_type}/{sample}/{type}",
    shell:
        """
        conda list &> {log};
        if [[ "{wildcards.data_type}" != "gencode" ]]; then
            salmon quant --ont -t {input.sirv_transcriptome} -a {input.reads} \
                    -o {params.outdir} -p {threads} -l A \
                    &>> {log}
        else
            salmon quant --ont -t {input.gencode_transcriptome} -a {input.reads} \
                    -o {params.outdir} -p {threads} -l A \
                    &>> {log}
        fi
        """


def get_downsampled_counts_paths(wildcards):
    method_to_path_suffix = {
        "run_bambu_lr": "counts_transcript.txt",
        "run_oarfish_lr": "",
        "run_isosceles_lr": "counts.tsv",
        "run_salmon_illumina": "quant.sf",
        "run_salmon_illumina_tpm": "quant.sf",
        "run_salmon_illumina_corrected": "quant.sf",
        "run_kallisto_long_lr": "abundance_1.tsv",
        "run_salmon_lr": "quant.sf",
        "run_isoquant_lr": "OUT/OUT.transcript_counts.tsv",
        "run_miniquant_lr": "abundance.tsv",
    }

    suffix = method_to_path_suffix[wildcards.method]

    if wildcards.method == "run_salmon_illumina_tpm":
        base_pattern = "results/quantify_bulk_downsampled/run_salmon_illumina/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"
    else:
        base_pattern = "results/quantify_bulk_downsampled/{method}/{subsample_number}_{read_number}/{data_type}/{sample}/{type}"

    if wildcards.method == "run_oarfish_lr":
        full_pattern = f"{base_pattern}/{{sample}}.quant"
    elif suffix:
        full_pattern = f"{base_pattern}/{suffix}"
    else:
        full_pattern = base_pattern

    return expand(
        full_pattern,
        sample=config["sample_names"],
        allow_missing=True,
        **wildcards,
    )


rule quantify_bulk_downsampled_format:
    input:
        counts=get_downsampled_counts_paths,
        gtf=lambda wildcards: (
            "results/setup/standardize_gtf_files/sirv_set_four.gtf"
            if wildcards.data_type != "gencode"
            else "results/setup/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf"
        ),
    output:
        out_transcript="results/quantify_bulk_downsampled/format/{method}/{subsample_number}_{read_number}/{type}/{data_type}/transcript_counts_formatted.tsv",
        out_gene="results/quantify_bulk_downsampled/format/{method}/{subsample_number}_{read_number}/{type}/{data_type}/gene_counts_formatted.tsv",
    log:
        "logs/quantify_bulk_downsampled/format/{method}/{subsample_number}_{read_number}/{type}/{data_type}/out.log",
    conda:
        "../envs/py/aggregate.yaml"
    threads: 1
    params:
        bambu_transcript_id_col_ix=1,
        bambu_tx_count_id_col_ix=3,
        isosceles_transcript_id_col_ix=1,
        isosceles_tx_count_id_col_ix=2,
        kallisto_transcript_id_col_ix=1,
        kallisto_tx_count_id_col_ix=4,
        salmon_transcript_id_col_ix=1,
        salmon_tx_count_id_col_ix=5,
        salmon_tpm_id_col_ix=4,
        sample_names=config["sample_names"],
    shell:
        """
        conda list &> {log};

        if [[ "{wildcards.method}" == "run_bambu_lr" ]] || [[ "{wildcards.method}" == "run_miniquant_lr" ]] || [[ "{wildcards.method}" == "run_oarfish_lr" ]]; then
            TRANSCRIPT_ID_COL={params.bambu_transcript_id_col_ix}
            COUNT_COL={params.bambu_tx_count_id_col_ix}
        elif [[ "{wildcards.method}" == "run_isosceles_lr" ]] || [[ "{wildcards.method}" == "run_isoquant_lr" ]]; then
            TRANSCRIPT_ID_COL={params.isosceles_transcript_id_col_ix}
            COUNT_COL={params.isosceles_tx_count_id_col_ix}
        elif [[ "{wildcards.method}" == "run_kallisto_long_lr" ]]; then
            TRANSCRIPT_ID_COL={params.kallisto_transcript_id_col_ix}
            COUNT_COL={params.kallisto_tx_count_id_col_ix}
        elif [[ "{wildcards.method}" == "run_salmon_lr" ]] || [[ "{wildcards.method}" == "run_salmon_illumina" ]]; then
            TRANSCRIPT_ID_COL={params.salmon_transcript_id_col_ix}
            COUNT_COL={params.salmon_tx_count_id_col_ix}
        elif [[ "{wildcards.method}" == "run_salmon_lr" ]] || [[ "{wildcards.method}" == "run_salmon_illumina_tpm" ]]; then
            TRANSCRIPT_ID_COL={params.salmon_transcript_id_col_ix}
            COUNT_COL={params.salmon_tpm_id_col_ix}
        elif [[ "{wildcards.method}" == "run_salmon_illumina_corrected" ]]; then
            TRANSCRIPT_ID_COL={params.isosceles_transcript_id_col_ix}
            COUNT_COL={params.isosceles_tx_count_id_col_ix}
        fi

        python workflow/scripts/py/aggregate_counts_debug.py \
                    --sample_paths {input.counts} \
                    --sample_names {params.sample_names} \
                    --gtf_annotation_path {input.gtf} \
                    --output_path_transcript {output.out_transcript} \
                    --output_path_gene {output.out_gene} \
                    --transcript_id_col_ix $TRANSCRIPT_ID_COL \
                    --count_col_ix $COUNT_COL &>> {log}
        """


rule quantify_bulk_joint_format:
    input:
        counts=expand(
            "results/quantify_bulk_joint/{{tech}}/run_miniquant_lr/{sample}/abundance.tsv",
            sample=config["sample_names"],
        ),
        gtf="results/setup/adjust_transcriptome_assembly_names/gencode.v45.primary_assembly.annotation.named.gtf",
    output:
        out_transcript="results/quantify_bulk_joint_format/format/{tech}/transcript_counts_formatted.tsv",
        out_gene="results/quantify_bulk_joint_format/format/{tech}/gene_counts_formatted.tsv",
    log:
        "logs/quantify_bulk_joint_format/{tech}/out.log",
    conda:
        "../envs/py/aggregate.yaml"
    threads: 1
    params:
        miniquant_transcript_id_col_ix=1,
        miniquant_tx_count_id_col_ix=3,
        sample_names=config["sample_names"],
    shell:
        """
        conda list &> {log};
        TRANSCRIPT_ID_COL={params.miniquant_transcript_id_col_ix}
        COUNT_COL={params.miniquant_tx_count_id_col_ix}
        python workflow/scripts/py/aggregate_counts_debug.py \
                    --sample_paths {input.counts} \
                    --sample_names {params.sample_names} \
                    --gtf_annotation_path {input.gtf} \
                    --output_path_transcript {output.out_transcript} \
                    --output_path_gene {output.out_gene} \
                    --transcript_id_col_ix $TRANSCRIPT_ID_COL \
                    --count_col_ix $COUNT_COL &>> {log}
        """
