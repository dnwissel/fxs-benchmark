configfile: "config/config.yaml"


rule download_genome:
    output:
        "results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    log:
        "logs/setup/download_genome/out.log",
    conda:
        "../envs/standalone/curl.yaml"
    params:
        url=config["genome_url"],
        download_path="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    shell:
        """
        conda list &> {log};
        wget -O {params.download_path} {params.url} &>> {log};
        gunzip {params.download_path}
        """


rule download_transcriptome:
    output:
        "results/setup/download_transcriptome/gencode.v49.primary_assembly.annotation.gtf.gz",
    log:
        "logs/setup/download_transcriptome/out.log",
    params:
        url=config["transcriptome_url"],
    shell:
        """
        wget -O {output} {params.url} &> {log}
        """


rule download_sirvome:
    output:
        sirv_set_four=directory(
            "results/setup/download_sirvome/SIRV_Set4_Norm_Sequences_20210507"
        ),
    log:
        "logs/setup/download_sirvome/out.log",
    conda:
        "../envs/standalone/curl.yaml"
    params:
        set_four_url=config["sirv_set_four_url"],
        set_four_download_object="results/setup/download_sirvome/set_four.zip",
        download_path="results/setup/download_sirvome/",
    shell:
        """
        conda list &> {log};
        curl -L {params.set_four_url} \
            --output {params.set_four_download_object} &>> {log};
        unzip -o {params.set_four_download_object} \
            -d {params.download_path} &>> {log}
        """


rule adjust_transcriptome_assembly_names:
    input:
        genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        transcriptome="results/setup/download_transcriptome/gencode.v49.primary_assembly.annotation.gtf.gz",
    output:
        "results/setup/adjust_transcriptome_assembly_names/gencode.v49.primary_assembly.annotation.named.gtf",
    log:
        "logs/setup/adjust_transcriptome_assembly_names/out.log",
    conda:
        "../envs/r/base.yaml"
    script:
        "../scripts/r/adjust_transcriptome_assembly_names.R"


rule adjust_sirv_names:
    input:
        "results/setup/download_sirvome/SIRV_Set4_Norm_Sequences_20210507",
    output:
        output_path_transcriptome="results/setup/adjust_sirv_names/sirv_set_four.gtf",
        output_path_genome="results/setup/adjust_sirv_names/sirv_set_four.fa",
    log:
        "logs/setup/adjust_sirv_names/out.log",
    conda:
        "../envs/r/base.yaml"
    params:
        transcriptome_name=config["sirv_transcriptome_name"],
        genome_name=config["sirv_genome_name"],
    script:
        "../scripts/r/adjust_sirv_names.R"


rule standardize_gtf_files:
    input:
        sirv_four_transcriptome="results/setup/adjust_sirv_names/sirv_set_four.gtf",
        gencode_transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v49.primary_assembly.annotation.named.gtf",
    output:
        sirv_four_transcriptome="results/setup/standardize_gtf_files/sirv_set_four.gtf",
        sirv_gff="results/setup/standardize_gtf_files/sirv_set_four.gff",
        gencode_transcriptome="results/setup/standardize_gtf_files/gencode.v49.primary_assembly.annotation.named.gtf",
        gencode_gff="results/setup/standardize_gtf_files/gencode.v49.primary_assembly.annotation.named.gff",
        gencode_transcriptome_gmap="results/setup/standardize_gtf_files/gencode_map.txt",
        gencode_transcriptome_gmap_headered="results/setup/standardize_gtf_files/gencode_map_headered.txt",
        sirv_four_transcriptome_gmap="results/setup/standardize_gtf_files/sirv_set_four_map.txt",
        sirv_four_transcriptome_gmap_headered="results/setup/standardize_gtf_files/sirv_set_four_map_headered.txt",
    log:
        "logs/setup/standardize_gtf_files/out.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -E {input.sirv_four_transcriptome} \
            -T -o {output.sirv_four_transcriptome} &>> {log};
        gffread -E {input.gencode_transcriptome} \
            -T -o {output.gencode_transcriptome} &>> {log};
        gffread {output.sirv_four_transcriptome} \
            --table transcript_id,gene_id \
            > {output.sirv_four_transcriptome_gmap} 2>> {log};
        gffread {output.gencode_transcriptome} \
            --table transcript_id,gene_id \
            > {output.gencode_transcriptome_gmap} 2>> {log};
        echo -e "transcript\tgene" | \
            cat - {output.gencode_transcriptome_gmap} \
                > {output.gencode_transcriptome_gmap_headered} \
                2>> {log};
        echo -e "transcript\tgene" | \
            cat - {output.sirv_four_transcriptome_gmap} \
                > {output.sirv_four_transcriptome_gmap_headered} \
                2>> {log};
        gffread -o {output.gencode_gff} \
            {output.gencode_transcriptome} &>> {log};
        gffread -o {output.sirv_gff} \
            {output.sirv_four_transcriptome} &>> {log}
        """


rule make_db_files_sirv:
    input:
        input_gtf="results/setup/standardize_gtf_files/sirv_set_four.gtf",
    output:
        output_db="results/setup/make_db_files/sirv_set_four.db",
    log:
        "logs/setup/make_db_files/sirv.log",
    conda:
        "../envs/py/gffutils.yaml"
    params:
        checklines=config["gffutils_checklines"],
    script:
        "../scripts/py/create_db_files.py"


rule make_db_files_gencode:
    input:
        input_gtf="results/setup/standardize_gtf_files/gencode.v49.primary_assembly.annotation.named.gtf",
    output:
        output_db="results/setup/make_db_files/gencode.v49.primary_assembly.annotation.named.db",
    log:
        "logs/setup/make_db_files/gencode.log",
    conda:
        "../envs/standalone/gffutils.yaml"
    params:
        checklines=config["gffutils_checklines"],
    script:
        "../scripts/py/create_db_files.py"


rule concatenate_genomes:
    input:
        sirv_four_genome="results/setup/adjust_sirv_names/sirv_set_four.fa",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        overall_genome="results/setup/concatenate_genomes/genome.fa",
    log:
        "logs/setup/concatenate_genomes/out.log",
    conda:
        "../envs/standalone/samtools.yaml"
    shell:
        """
        conda list &> {log};
        cat {input.gencode_genome} {input.sirv_four_genome} > \
            {output.overall_genome} 2>> {log};
        samtools faidx {output.overall_genome} &>> {log}
        """


rule concatenate_transcriptomes:
    input:
        sirv_four_transcriptome="results/setup/standardize_gtf_files/sirv_set_four.gtf",
        gencode_transcriptome="results/setup/standardize_gtf_files/gencode.v49.primary_assembly.annotation.named.gtf",
    output:
        overall_transcriptome="results/setup/concatenate_transcriptomes/transcriptome.gtf",
    log:
        "logs/setup/concatenate_transcriptomes/out.log",
    shell:
        """
        cat {input.gencode_transcriptome} {input.sirv_four_transcriptome} \
            > {output.overall_transcriptome} 2> {log}
        """


rule convert_gtfs_to_beds_gencode:
    input:
        gencode_gtf="results/setup/standardize_gtf_files/gencode.v49.primary_assembly.annotation.named.gtf",
        sirv_gtf="results/setup/standardize_gtf_files/sirv_set_four.gtf",
        overall_gtf="results/setup/concatenate_transcriptomes/transcriptome.gtf",
    output:
        gencode_bed="results/setup/convert_gtfs_to_beds/gencode.bed",
        sirv_bed="results/setup/convert_gtfs_to_beds/sirv.bed",
        overall_bed="results/setup/convert_gtfs_to_beds/overall.bed",
    log:
        "logs/setup/convert_gtfs_to_beds/out.log",
    conda:
        "../envs/standalone/minimap2.yaml"
    shell:
        """
        conda list &> {log};
        paftools.js gff2bed {input.gencode_gtf} > \
            {output.gencode_bed} 2>> {log};
        paftools.js gff2bed {input.sirv_gtf} > \
            {output.sirv_bed} 2>> {log};
        paftools.js gff2bed {input.overall_gtf} > \
            {output.overall_bed} 2>> {log};
        """


rule extract_transcriptomes:
    input:
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v49.primary_assembly.annotation.named.gtf",
        joint_genome="results/setup/concatenate_genomes/genome.fa",
        joint_transcriptome="results/setup/concatenate_transcriptomes/transcriptome.gtf",
        sirv_genome="results/setup/adjust_sirv_names/sirv_set_four.fa",
        sirv_transcriptome="results/setup/adjust_sirv_names/sirv_set_four.gtf",
    output:
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
        joint_transcriptome="results/setup/extract_transcriptomes/joint_transcriptome.fa",
        sirv_transcriptome="results/setup/extract_transcriptomes/sirv_transcriptome.fa",
    log:
        "logs/setup/extract_transcriptomes/out.log",
    conda:
        "../envs/standalone/gffread.yaml"
    shell:
        """
        conda list &> {log};
        gffread -w {output.gencode_transcriptome} \
            -g {input.gencode_genome} \
            {input.gencode_transcriptome} &>> {log};
        gffread -w {output.sirv_transcriptome} \
            -g {input.sirv_genome} \
            {input.sirv_transcriptome} &>> {log};
        gffread -w {output.joint_transcriptome} \
            -g {input.joint_genome} \
            {input.joint_transcriptome} &>> {log}
        """


rule index_star:
    input:
        genome="results/setup/concatenate_genomes/genome.fa",
        transcriptome="results/setup/concatenate_transcriptomes/transcriptome.gtf",
    output:
        directory("results/setup/index_star"),
    log:
        "logs/setup/index_star/out.log",
    conda:
        "../envs/standalone/star.yaml"
    threads: config["align_map_bam_threads"]
    params:
        read_length=config["illumina_read_length"],
    shell:
        """
        conda list &> {log};
        STAR --runMode genomeGenerate --runThreadN {threads} \
            --genomeDir {output} --genomeFastaFiles {input.genome} \
            --sjdbGTFfile {input.transcriptome} \
            --sjdbOverhang {params.read_length} &>> {log}
        """


rule create_decoys_salmon_index:
    input:
        sirv_genome="results/setup/adjust_sirv_names/sirv_set_four.fa",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        gencode_transcriptome="results/setup/extract_transcriptomes/gencode_transcriptome.fa",
        sirv_transcriptome="results/setup/extract_transcriptomes/sirv_transcriptome.fa",
    output:
        sirv_decoys="results/setup/create_decoys_salmon_index/sirv_decoys.txt",
        gencode_decoys="results/setup/create_decoys_salmon_index/gencode_decoys.txt",
        gencode_decoyed_transcriptome="results/setup/create_decoys_salmon_index/gencode_decoyed_sequence.fa",
        sirv_decoyed_transcriptome="results/setup/create_decoys_salmon_index/sirv_decoyed_sequence.fa",
    log:
        "logs/setup/create_decoys_salmon_index/out.log",
    params:
        outdir="results/setup/create_decoys_salmon_index",
    shell:
        """
        touch {log};
        mkdir -p {params.outdir} &>> {log};
        grep "^>" {input.sirv_genome} | cut -d " " -f 1 \
            > {output.sirv_decoys} 2>> {log}
        grep "^>" {input.gencode_genome} | cut -d " " -f 1 \
            > {output.gencode_decoys} 2>> {log};
        sed -i.bak -e 's/>//g' {output.sirv_decoys} 2>> {log};
        sed -i.bak -e 's/>//g' {output.gencode_decoys} 2>> {log};
        cat {input.gencode_transcriptome} {input.gencode_genome} \
            > {output.gencode_decoyed_transcriptome} 2>> {log};
        cat {input.sirv_transcriptome} {input.sirv_genome} \
            > {output.sirv_decoyed_transcriptome} 2>> {log}
        """


rule create_salmon_index:
    input:
        decoys="results/setup/create_decoys_salmon_index/{type}_decoys.txt",
        decoyed_transcriptome="results/setup/create_decoys_salmon_index/{type}_decoyed_sequence.fa",
    output:
        directory("results/setup/create_salmon_index/{type}"),
    log:
        "logs/setup/create_salmon_index/{type}.log",
    conda:
        "../envs/standalone/salmon.yaml"
    threads: config["align_map_bam_threads"]
    params:
        k=config["salmon_index_k"],
    shell:
        """
        conda list &> {log};
        salmon index -t {input.decoyed_transcriptome} \
            -i {output} -d {input.decoys} --keepDuplicates \
            -k {params.k} -p {threads}  &>> {log}
        """


rule compile_lr_kallisto:
    output:
        directory("results/setup/compile_lr_kallisto"),
    log:
        "/home/dwissel/projects/fxs-benchmark/logs/setup/compile_lr_kallisto/out.log",
    conda:
        "../envs/standalone/cxx.yaml"
    shell:
        """
        conda list &> {log};
        mkdir {output}; 
        cd {output};
        git clone https://github.com/pachterlab/kallisto &>> {log};
        cd kallisto &>> {log};
        mkdir build;
        cd build &> {log};
        cmake .. -DMAX_KMER_SIZE=64 &> {log};
        make &>> {log};
        """


rule create_lr_kallisto_index:
    input:
        kallisto_binary="results/setup/compile_lr_kallisto",
        transcriptome="results/setup/extract_transcriptomes/{type}_transcriptome.fa",
    output:
        f"results/setup/create_lr_kallisto_index/{{type}}_k-{config['lr_kallisto_index_k']}.idx",
    log:
        "logs/setup/create_lr_kallisto_index/{type}.log",
    threads: config["align_map_bam_threads"]
    params:
        k=config["lr_kallisto_index_k"],
    shell:
        """
        touch {log};
        {input.kallisto_binary}/kallisto/build/src/kallisto index \
            -k {params.k} -t {threads} -i {output} \
            {input.transcriptome} &>> {log}
        """


rule create_simpleaf_index:
    input:
        genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
        transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v49.primary_assembly.annotation.named.gtf",
    output:
        directory("results/setup/create_simpleaf_index"),
        "results/setup/create_simpleaf_index/index/piscem_idx.sshash",
        "results/setup/create_simpleaf_index/index/piscem_idx.refinfo",
        "results/setup/create_simpleaf_index/index/t2g_3col.tsv",
    log:
        "logs/setup/create_simpleaf_index/out.log",
    conda:
        "../envs/standalone/simpleaf.yaml"
    threads: config["align_map_bam_threads"]
    params:
        outdir="results/setup/create_simpleaf_index",
    shell:
        """
        conda list &> {log};
        export ALEVIN_FRY_HOME="{output}";
        simpleaf set-paths &>> {log};
        ulimit -n 2048 &>> {log};
        simpleaf index \
            --output {params.outdir} \
            --fasta {input.genome} \
            --gtf {input.transcriptome} \
            --threads {threads} \
            --work-dir . &>> {log}
        """


rule install_bambu:
    output:
        "results/setup/install_bambu/done.txt",
    log:
        "logs/setup/install_bambu/out.log",
    conda:
        "../envs/r/bambu.yaml"
    params:
        version=config["bambu_version"],
    script:
        "../scripts/r/install_bambu.R"


rule setup_install_isosceles:
    output:
        "results/setup/install_isosceles/done.txt",
    log:
        "logs/setup/install_isosceles/out.log",
    conda:
        "../envs/r/isosceles.yaml"
    params:
        version=config["isosceles_version"],
    script:
        "../scripts/r/install_isosceles.R"


rule filter_tusco:
    input:
        transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v49.primary_assembly.annotation.named.gtf",
        tusco="config/tusco_human.tsv",
    output:
        transcriptome_without_tusco="results/setup/filter_tusco/gencode_no_tusco.gtf",
        tusco="results/setup/filter_tusco/tusco.gtf",
    log:
        "logs/setup/filter_tusco/out.log",
    conda:
        "../envs/r/r_base.yaml"
    script:
        "../scripts/r/filter_tusco.R"


rule filter_gencode:
    input:
        transcriptome="results/setup/adjust_transcriptome_assembly_names/gencode.v49.primary_assembly.annotation.named.gtf",
        illumina="results/quantify_bulk_downsampled/format/run_salmon_illumina_corrected/1_15000000.0/illumina/gencode/transcript_counts_formatted.tsv",
        kinnex="results/quantify_bulk_downsampled/format/run_oarfish_lr/1_15000000.0/pb/gencode/transcript_counts_formatted.tsv",
        ont="results/quantify_bulk_downsampled/format/run_oarfish_lr/1_15000000.0/ont/gencode/transcript_counts_formatted.tsv",
    output:
        "results/setup/filter_gencode/gencode_removed.gtf",
        "results/setup/filter_gencode/gencode_removed_only.gtf",
    log:
        "logs/setup/filter_gencode/out.log",
    conda:
        "../envs/r/base_edger.yaml"
    script:
        "../scripts/r/filter_gencode.R"
