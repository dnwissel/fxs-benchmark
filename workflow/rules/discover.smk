rule discover_sc_run_bambu:
    input:
        "results/setup/install_bambu/done.txt",
        "results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam.bai",
        reads="results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam",
        gencode_transcriptome="results/setup/filter_tusco/gencode_no_tusco.gtf",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/discover_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
    log:
        "logs/discover_sc_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/r/bambu.yaml"
    threads: config["quantify_threads"]
    script:
        "../scripts/r/run_bambu_discovery.R"


rule discover_bulk_run_bambu:
    input:
        "results/setup/install_bambu/done.txt",
        reads="results/quantify_downsampled/align_genome_bulk_lr/{subsample_number}_{number_to_sample}/gencode/{sample}/{tech}/{sample}.aligned.sorted.bam",
        gencode_transcriptome="results/setup/filter_tusco/gencode_no_tusco.gtf",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/discover_bulk_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
    log:
        "logs/discover_bulk_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/r/bambu.yaml"
    threads: config["quantify_threads"]
    script:
        "../scripts/r/run_bambu_discovery.R"


rule discover_sc_run_bambu_real:
    input:
        "results/setup/install_bambu/done.txt",
        "results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam.bai",
        reads="results/preprocess/adjust_lr_sc_genome/{tech}/{subsample_number}_{number_to_sample}/{sample}.bam",
        gencode_transcriptome="results/setup/filter_gencode/gencode_removed.gtf",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/discover_sc_run_bambu_real/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
    log:
        "logs/discover_sc_run_bambu_real/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/r/bambu.yaml"
    threads: config["quantify_threads"]
    script:
        "../scripts/r/run_bambu_discovery.R"


rule discover_bulk_run_bambu_real:
    input:
        "results/setup/install_bambu/done.txt",
        reads="results/quantify_downsampled/align_genome_bulk_lr/{subsample_number}_{number_to_sample}/gencode/{sample}/{tech}/{sample}.aligned.sorted.bam",
        gencode_transcriptome="results/setup/filter_gencode/gencode_removed.gtf",
        gencode_genome="results/setup/download_genome/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna",
    output:
        "results/discover_bulk_run_bambu_real/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
    log:
        "logs/discover_bulk_run_bambu_real/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/r/bambu.yaml"
    threads: config["quantify_threads"]
    script:
        "../scripts/r/run_bambu_discovery.R"


rule discover_filter_bambu:
    input:
        transcriptome="results/discover_{level}_run_bambu/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
        tusco="results/setup/filter_tusco/tusco.gtf",
    output:
        "results/discover_filter_bambu/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
    log:
        "logs/discover_filter_bambu/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/r/r_base.yaml"
    threads: config["quantify_threads"]
    script:
        "../scripts/r/filter_bambu.R"


rule discover_filter_bambu_novel:
    input:
        transcriptome="results/discover_{level}_run_bambu_real/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
        tusco="results/setup/filter_tusco/tusco.gtf",
    output:
        "results/discover_filter_bambu_novel/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
    log:
        "logs/discover_filter_bambu_novel/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/r/r_base.yaml"
    threads: config["quantify_threads"]
    script:
        "../scripts/r/filter_bambu_novel.R"


rule evaluate_performance_tusco:
    input:
        transcriptome="results/discover_filter_bambu/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
        ground_truth="results/setup/filter_tusco/tusco.gtf",
    output:
        "results/evaluate_performance_tusco/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.stats",
    log:
        "logs/evaluate_performance_tusco/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/gffcompare.yaml"
    params:
        output="results/evaluate_performance_tusco/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}",
    shell:
        """
        conda list > {log};
        gffcompare -r {input.ground_truth} -o {params.output} \
            {input.transcriptome} &>> {log}
        """


rule evaluate_performance_real:
    input:
        transcriptome="results/discover_filter_bambu_novel/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.gtf",
        ground_truth="results/setup/filter_gencode/gencode_removed_only.gtf",
    output:
        "results/evaluate_performance_real/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.stats",
    log:
        "logs/evaluate_performance_real/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}.log",
    conda:
        "../envs/standalone/gffcompare.yaml"
    params:
        output="results/evaluate_performance_real/{level}/{tech}/{subsample_number}_{number_to_sample}/{sample}",
    shell:
        """
        conda list > {log};
        gffcompare -r {input.ground_truth} -o {params.output} \
            {input.transcriptome} &>> {log}
        """
