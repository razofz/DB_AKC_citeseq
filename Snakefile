import os


configfile: "config.yaml"


rule produce_unprocessed_seurat_object:
    input:
        filtered_h5=config["raw_dir"] + "filtered_feature_bc_matrix.h5"
    output:
        seurat_object_unprocessed=config["interim_dir"] + "seurat_object_unprocessed.rds",
        cell_whitelist=config["interim_dir"] + "cell_whitelist.tsv",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/demultiplexing/pre-cellhashr-seuratobject.R"


raw_counts_dir = config["interim_dir"] + "raw_counts_htos/"


rule write_raw_counts_hto:
    input:
        raw_counts=config["raw_dir"] + "raw_feature_bc_matrix/"
    output:
        barcodes=raw_counts_dir + "barcodes.tsv",
        features=raw_counts_dir + "features.tsv",
        mtx=raw_counts_dir + "matrix.mtx",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/demultiplexing/pre-cellhashr.R"


rule gzip_raw_counts_hto:
    input:
        barcodes=rules.write_raw_counts_hto.output.barcodes,
        features=rules.write_raw_counts_hto.output.features,
        mtx=rules.write_raw_counts_hto.output.mtx,
    output:
        barcodes=raw_counts_dir + "barcodes.tsv.gz",
        features=raw_counts_dir + "features.tsv.gz",
        mtx=raw_counts_dir + "matrix.mtx.gz",
    conda:
        "envs/DB_AKC_gzip.yaml"
    shell:
        "for file in {input}; do gzip --keep $file; done"


cellhashr_dir = config["processed_dir"] + "cellhashr_output/"


rule run_cellhashr:
    input:
        barcodes=rules.gzip_raw_counts_hto.output.barcodes,
        features=rules.gzip_raw_counts_hto.output.features,
        mtx=rules.gzip_raw_counts_hto.output.mtx,
        cell_whitelist=rules.produce_unprocessed_seurat_object.output.cell_whitelist,
    output:
        metrics_file_ProcessCountMatrix=cellhashr_dir + "metrics_file_ProcessCountMatrix.tsv",
        metrics_file_GenerateCellHashingCalls=cellhashr_dir + "metrics_file_GenerateCellHashingCall.tsv",
        rplots_pdf=cellhashr_dir + "hashing_plots.pdf",
        demultiplexing_results=cellhashr_dir + "demultiplexing_results.tsv",
    container:
        "docker://ghcr.io/bimberlab/cellhashr:latest"
    resources:
        mem_mb=config["max_memory_to_use_in_gb"] * 1000,
    script:
        "src/smk/demultiplexing/run-cellhashr.R"


HTO_out_dir = config["processed_dir"] + "HTO_demultiplexed/"


rule plot_hto_distribution:
    input:
        seurat_object_unprocessed=rules.produce_unprocessed_seurat_object.output.seurat_object_unprocessed,
        demultiplexing_results=rules.run_cellhashr.output.demultiplexing_results,
    output:
        hto_distribution_cellhashr=HTO_out_dir + "HTO_distribution.svg",
        seurat_object_demultiplexed=HTO_out_dir + "seurat_object_demultiplexed.rds",
    script:
        "src/smk/demultiplexing/post-cellhashr.R"


