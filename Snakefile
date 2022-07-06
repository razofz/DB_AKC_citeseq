import os


configfile: "config.yaml"


rule produce_unprocessed_seurat_object:
    input:
        filtered_h5=config["raw_dir"] + "filtered_feature_bc_matrix.h5",
    output:
        seurat_object_unprocessed=config["interim_dir"]
        + "seurat_object_unprocessed.rds",
        cell_whitelist=config["interim_dir"] + "cell_whitelist.tsv",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/demultiplexing/pre_cellhashr_seuratobject.R"


raw_counts_dir = config["interim_dir"] + "raw_counts_htos/"


rule write_raw_counts_hto:
    input:
        raw_counts=config["raw_dir"] + "raw_feature_bc_matrix/",
    output:
        barcodes=raw_counts_dir + "barcodes.tsv",
        features=raw_counts_dir + "features.tsv",
        mtx=raw_counts_dir + "matrix.mtx",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/demultiplexing/pre_cellhashr.R"


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
        metrics_file_ProcessCountMatrix=cellhashr_dir
        + "metrics_file_ProcessCountMatrix.tsv",
        metrics_file_GenerateCellHashingCalls=cellhashr_dir
        + "metrics_file_GenerateCellHashingCall.tsv",
        rplots_pdf=cellhashr_dir + "hashing_plots.pdf",
        demultiplexing_results=cellhashr_dir + "demultiplexing_results.tsv",
    container:
        "docker://ghcr.io/bimberlab/cellhashr:latest"
    resources:
        mem_mb=config["max_memory_to_use_in_gb"] * 1000,
    script:
        "src/smk/demultiplexing/run_cellhashr.R"


HTO_out_dir = config["processed_dir"] + "HTO_demultiplexed/"


rule post_cellhashr:
    input:
        seurat_object_unprocessed=rules.produce_unprocessed_seurat_object.output.seurat_object_unprocessed,
        demultiplexing_results=rules.run_cellhashr.output.demultiplexing_results,
    output:
        seurat_object_demultiplexed=HTO_out_dir + "seurat_object_demultiplexed.rds",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/demultiplexing/post_cellhashr.R"


rule plot_hto_distribution:
    input:
        seurat_object=rules.plot_hto_distribution.output.seurat_object_demultiplexed,
    output:
        hto_distribution_cellhashr=report(
            HTO_out_dir + "HTO_distribution.svg",
            category="Demultiplexing"
        ),
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/plot_hto_distribution.R"


filtering_dir = config["processed_dir"] + "seurat_filtering/"


rule plot_prefiltering:
    input:
        seurat_object=rules.plot_hto_distribution.output.seurat_object_demultiplexed,
    output:
        RNA_vln=report(filtering_dir + "pre_RNA.svg", category="Filtering"),
        RNA_vln_log=report(filtering_dir + "pre_RNA_log.svg", category="Filtering"),
        ADT_vln=report(filtering_dir + "pre_ADT.svg", category="Filtering"),
        ADT_vln_log=report(filtering_dir + "pre_ADT_log.svg", category="Filtering"),
        HTO_vln=report(filtering_dir + "pre_HTO.svg", category="Filtering"),
        HTO_vln_log=report(filtering_dir + "pre_HTO_log.svg", category="Filtering"),
        mito_vln=report(filtering_dir + "pre_mito.svg", category="Filtering"),
        mito_vln_log=report(filtering_dir + "pre_mito_log.svg", category="Filtering"),
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/plot_prefiltering.R"


rule filtering:
    input:
        seurat_object=rules.post_cellhashr.output.seurat_object,
    output:
        seurat_object=filtering_dir + "seurat_object_filtered.rds",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/seurat_filtering.R"


rule filtering:
    input:
        seurat_object=rules.filtering.output.seurat_object,
    output:
        RNA_vln=report(filtering_dir + "post_RNA.svg", category="Filtering"),
        RNA_vln_log=report(filtering_dir + "post_RNA_log.svg", category="Filtering"),
        ADT_vln=report(filtering_dir + "post_ADT.svg", category="Filtering"),
        ADT_vln_log=report(filtering_dir + "post_ADT_log.svg", category="Filtering"),
        HTO_vln=report(filtering_dir + "post_HTO.svg", category="Filtering"),
        HTO_vln_log=report(filtering_dir + "post_HTO_log.svg", category="Filtering"),
        mito_vln=report(filtering_dir + "post_mito.svg", category="Filtering"),
        mito_vln_log=report(filtering_dir + "post_mito_log.svg", category="Filtering"),
        hto_distribution_cellhashr=report(
            HTO_out_dir + "HTO_distribution_post_filtering.svg",
            category="Demultiplexing",
        ),
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/plot_post_filtering.R"


seurat_processing_all_dir = config["processed_dir"] + "seurat_processing_all/"


rule seurat_processing_all:
    input:
        seurat_object=rules.filtering.output.seurat_object,
    output:
        hvgs=seurat_processing_all_dir + "highly_variable_genes.tsv",
        seurat_object=seurat_processing_all_dir + "seurat_object.rds",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/seurat_processing_all.R"


rule seurat_umap_plotting:
    input:
        seurat_object=rules.seurat_processing_all.output.seurat_object,
    output:
        umap_clusters=report(
            seurat_processing_all_dir + "umap_clusters.svg",
            category="Seurat processing"
        ),
        umap_clusters_unlabelled=report(
            seurat_processing_all_dir + "umap_clusters_unlabelled.svg",
            category="Seurat processing",
        ),
        umap_hto=report(
            seurat_processing_all_dir + "umap_hto.svg",
            category="Seurat processing"
        ),
        umap_sample=report(
            seurat_processing_all_dir + "umap_sample.svg",
            category="Seurat processing"
        ),
        umap_incubation_method=report(
            seurat_processing_all_dir + "umap_incubation_method.svg",
            category="Seurat processing",
        ),
        umap_buffer_treatment=report(
            seurat_processing_all_dir + "umap_buffer_treatment.svg",
            category="Seurat processing",
        ),
        umap_cellcycle=report(
            seurat_processing_all_dir + "umap_cellcycle.svg",
            category="Seurat processing",
        ),
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/seurat_umap_plotting.R"


deg_dir = config["processed_dir"] + "deg/"


rule seurat_deg:
    input:
        seurat_object=rules.seurat_processing_all.output.seurat_object,
    output:
        deg_S_triptolides=deg_dir + "deg_S_triptolides.tsv",
        deg_G1_triptolides=deg_dir + "deg_G1_triptolides.tsv",
        deg_G2M_triptolides=deg_dir + "deg_G2M_triptolides.tsv",
        deg_S_no_triptolides=deg_dir + "deg_S_no_triptolides.tsv",
        deg_G1_no_triptolides=deg_dir + "deg_G1_no_triptolides.tsv",
        deg_G2M_no_triptolides=deg_dir + "deg_G2M_no_triptolides.tsv",
        deg_S_ice_vs_37c_t=deg_dir + "deg_S_ice_vs_37c_t.tsv",
        deg_G1_ice_vs_37c_t=deg_dir + "deg_G1_ice_vs_37c_t.tsv",
        deg_G2M_ice_vs_37c_t=deg_dir + "deg_G2M_ice_vs_37c_t.tsv",
        deg_S_ice_vs_37c_no_t=deg_dir + "deg_S_ice_vs_37c_no_t.tsv",
        deg_G1_ice_vs_37c_no_t=deg_dir + "deg_G1_ice_vs_37c_no_t.tsv",
        deg_G2M_ice_vs_37c_no_t=deg_dir + "deg_G2M_ice_vs_37c_no_t.tsv",
        deg_S_37c_no_t_vs_37c_t=deg_dir + "deg_S_37c_no_t_vs_37c_t.tsv",
        deg_G1_37c_no_t_vs_37c_t=deg_dir + "deg_G1_37c_no_t_vs_37c_t.tsv",
        deg_G2M_37c_no_t_vs_37c_t=deg_dir + "deg_G2M_37c_no_t_vs_37c_t.tsv",
        # seurat_object=deg_dir + "seurat_object_with_deg.rds",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/seurat_DEG.R"


seurat_processing_ice_dir = config["processed_dir"] + "seurat_processing_ice/"


rule ice_seurat_processing:
    input:
        seurat_object_all=rules.filtering.output.seurat_object,
    output:
        seurat_object_ice=seurat_processing_ice_dir + "seurat_object.rds",
    script:
        "src/smk/seurat_processing_ice.R"


