import os


configfile: "config.yaml"


raw_counts_dir = config["interim_dir"] + "raw_counts_htos/"
cellhashr_dir = config["processed_dir"] + "cellhashr_output/"
HTO_out_dir = config["processed_dir"] + "HTO_demultiplexed/"
filtering_dir = config["processed_dir"] + "seurat_filtering/"
seurat_processing_all_dir = config["processed_dir"] + "seurat_processing_all/"
deg_dir = config["processed_dir"] + "deg/"
zarr_dir = config["processed_dir"] + "zarr_files/"


rule gather:
    input:
        filtering_dir + "pre_RNA.svg",
        filtering_dir + "post_RNA.svg",
        seurat_processing_all_dir + "umap_cellcycle.svg",
        deg_dir + "deg_S_ice_t_vs_37c_t.tsv",
        config["processed_dir"] + "seurat_object_w_stress_sig.rds",
        config["processed_dir"] + "umap_sample.svg",
        seurat_processing_all_dir + "umap_embeddings.csv",


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


rule plot_prefiltering:
    input:
        seurat_object=rules.post_cellhashr.output.seurat_object_demultiplexed,
    output:
        RNA_vln=report(filtering_dir + "pre_RNA.svg", category="Filtering"),
        RNA_vln_log=report(filtering_dir + "pre_RNA_log.svg", category="Filtering"),
        ADT_vln=report(filtering_dir + "pre_ADT.svg", category="Filtering"),
        ADT_vln_log=report(filtering_dir + "pre_ADT_log.svg", category="Filtering"),
        HTO_vln=report(filtering_dir + "pre_HTO.svg", category="Filtering"),
        HTO_vln_log=report(filtering_dir + "pre_HTO_log.svg", category="Filtering"),
        mito_vln=report(filtering_dir + "pre_mito.svg", category="Filtering"),
        mito_vln_log=report(filtering_dir + "pre_mito_log.svg", category="Filtering"),
        hto_distribution_cellhashr=report(
            HTO_out_dir + "HTO_distribution.svg", category="Demultiplexing"
        ),
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/visualization/plot_prefiltering.R"


rule filtering:
    input:
        seurat_object=rules.post_cellhashr.output.seurat_object_demultiplexed,
    output:
        seurat_object=filtering_dir + "seurat_object_filtered.rds",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/seurat_filtering.R"


rule plot_filtered:
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
        "src/smk/visualization/plot_filtered.R"


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


rule plot_umap:
    input:
        seurat_object=rules.seurat_processing_all.output.seurat_object,
    output:
        umap_clusters=report(
            seurat_processing_all_dir + "umap_clusters.svg",
            category="Seurat processing",
        ),
        umap_clusters_unlabelled=report(
            seurat_processing_all_dir + "umap_clusters_unlabelled.svg",
            category="Seurat processing",
        ),
        umap_hto=report(
            seurat_processing_all_dir + "umap_hto.svg", category="Seurat processing"
        ),
        umap_sample=report(
            seurat_processing_all_dir + "umap_sample.svg", category="Seurat processing"
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
        "src/smk/visualization/plot_umap.R"


rule seurat_deg:
    input:
        seurat_object=rules.seurat_processing_all.output.seurat_object,
    output:
        deg_S_ice_t_vs_37c_t=deg_dir + "deg_S_ice_t_vs_37c_t.tsv",
        deg_G1_ice_t_vs_37c_t=deg_dir + "deg_G1_ice_t_vs_37c_t.tsv",
        deg_G2M_ice_t_vs_37c_t=deg_dir + "deg_G2M_ice_t_vs_37c_t.tsv",
        deg_S_ice_no_t_vs_37c_no_t=deg_dir + "deg_S_ice_no_t_vs_37c_no_t.tsv",
        deg_G1_ice_no_t_vs_37c_no_t=deg_dir + "deg_G1_ice_no_t_vs_37c_no_t.tsv",
        deg_G2M_ice_no_t_vs_37c_no_t=deg_dir + "deg_G2M_ice_no_t_vs_37c_no_t.tsv",
        deg_S_ice_vs_37c_t=deg_dir + "deg_S_ice_vs_37c_t.tsv",
        deg_G1_ice_vs_37c_t=deg_dir + "deg_G1_ice_vs_37c_t.tsv",
        deg_G2M_ice_vs_37c_t=deg_dir + "deg_G2M_ice_vs_37c_t.tsv",
        deg_S_ice_vs_37c_no_t=deg_dir + "deg_S_ice_vs_37c_no_t.tsv",
        deg_G1_ice_vs_37c_no_t=deg_dir + "deg_G1_ice_vs_37c_no_t.tsv",
        deg_G2M_ice_vs_37c_no_t=deg_dir + "deg_G2M_ice_vs_37c_no_t.tsv",
        deg_S_37c_no_t_vs_37c_t=deg_dir + "deg_S_37c_no_t_vs_37c_t.tsv",
        deg_G1_37c_no_t_vs_37c_t=deg_dir + "deg_G1_37c_no_t_vs_37c_t.tsv",
        deg_G2M_37c_no_t_vs_37c_t=deg_dir + "deg_G2M_37c_no_t_vs_37c_t.tsv",
        deg_S_ice_no_t_vs_ice_t=deg_dir + "deg_S_ice_no_t_vs_ice_t.tsv",
        deg_G1_ice_no_t_vs_ice_t=deg_dir + "deg_G1_ice_no_t_vs_ice_t.tsv",
        deg_G2M_ice_no_t_vs_ice_t=deg_dir + "deg_G2M_ice_no_t_vs_ice_t.tsv",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/seurat_DEG.R"


rule seurat_stress_signature:
    input:
        seurat_object=rules.seurat_processing_all.output.seurat_object,
        stress_signature=config["raw_base_dir"] + "stress_signature.tsv",
    output:
        seurat_object=config["processed_dir"] + "seurat_object_w_stress_sig.rds",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/add_gene_signatures.R"


rule plot_stress_signature:
    input:
        seurat_object=rules.seurat_stress_signature.output.seurat_object,
    output:
        umap_sample=config["processed_dir"] + "umap_sample.svg",
        umap_stress_sig_classification=config["processed_dir"] +
        "umap_stress_sig_classification.svg",
        # umap_stress_sig_score=config["processed_dir"] + "umap_stress_sig_score.svg",
        umap_stress_sig_score_q100=config["processed_dir"] +
        "umap_stress_sig_score_q100.svg",
        umap_stress_sig_score_q99=config["processed_dir"] +
        "umap_stress_sig_score_q99.svg",
        umap_stress_sig_score_blue=config["processed_dir"] +
        "umap_stress_sig_score_blue.svg",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/visualization/plot_stress_signature.R"


rule export_seurat_data:
    input:
        seurat_object=rules.seurat_processing_all.output.seurat_object,
    output:
        umap_embeddings=seurat_processing_all_dir + "umap_embeddings.csv",
        metadata=seurat_processing_all_dir + "metadata.csv",
    conda:
        "envs/DB_AKC_R.yaml"
    script:
        "src/smk/seurat_export_data.R"


# rule scarf_import_and_split_on_incubation:
#     input:
#         filtered_h5=config["raw_dir"] + "filtered_feature_bc_matrix.h5",
#         seurat_metadata=rules.export_seurat_data.output.metadata,
#     output:
#         zarr_all=directory(zarr_dir + "all.zarr"),
#         zarr_ice=directory(zarr_dir + "ice.zarr"),
#         zarr_37c=directory(zarr_dir + "37c.zarr"),
#     conda:
#         "envs/DB_AKC_scarf.yaml"
#     params:
#         threads=config["n_threads"],
#     script:
#         "src/smk/scarf_import_and_split_on_incubation.py"


# rule scarf_make_ice_reference:
#     input:
#         zarr_ice=rules.scarf_import_and_split_on_incubation.output.zarr_ice,
#     output:
#         ice_reference_marker=touch(".smk_markers/ice_reference_marker.done"),
#     conda:
#         "envs/DB_AKC_scarf.yaml"
#     params:
#         threads=config["n_threads"],
#     script:
#         "src/smk/scarf_make_ice_reference.py"


# scarf_dir = config["processed_dir"] + "scarf/"


# rule scarf_plot_ice_reference:
#     input:
#         zarr_ice=rules.scarf_import_and_split_on_incubation.output.zarr_ice,
#         ice_reference_marker=rules.scarf_make_ice_reference.output.ice_reference_marker,
#     output:
#         plot_umap_clusters=scarf_dir + "ice_umap_clusters.svg",
#         # plot_umap_buffer_treatment=scarf_dir + "ice_umap_buffer_treatment.svg",
#         # plot_umap_hto=scarf_dir + "ice_umap_hto.svg",
#         plot_umap_phase=scarf_dir + "ice_umap_phase.svg",
#         plot_umap_sample=scarf_dir + "ice_umap_sample.svg",
#     conda:
#         "envs/DB_AKC_scarf.yaml"
#     params:
#         threads=config["n_threads"],
#     script:
#         "src/smk/scarf_plot_ice_reference.py"


# rule scarf_make_37c_targets:
#     input:
#         zarr_37c=rules.scarf_import_and_split_on_incubation.output.zarr_37c,
#     output:
#         zarr_37c_t=directory(zarr_dir + "37c_t.zarr"),
#         zarr_37c_no_t=directory(zarr_dir + "37c_no_t.zarr"),
#         # 37_targets_marker=touch(".smk_markers/37_targets_marker.done"),
#     conda:
#         "envs/DB_AKC_scarf.yaml"
#     params:
#         threads=config["n_threads"],
#     script:
#         "src/smk/scarf_make_37c_targets.py"


# rule scarf_plot_37c_targets:
#     input:
#         zarr_37c_t=rules.scarf_make_37c_targets.output.zarr_37c_t,
#         zarr_37c_no_t=rules.scarf_make_37c_targets.output.zarr_37c_no_t,
#     output:
#         plot_umap_clusters_37c_t=scarf_dir + "umap_clusters_37c_t.svg",
#         plot_umap_phase_37c_t=scarf_dir + "umap_phase_37c_t.svg",
#         plot_umap_sample_37c_t=scarf_dir + "umap_sample_37c_t.svg",
#         plot_umap_clusters_37c_no_t=scarf_dir + "umap_clusters_37c_no_t.svg",
#         plot_umap_phase_37c_no_t=scarf_dir + "umap_phase_37c_no_t.svg",
#         plot_umap_sample_37c_no_t=scarf_dir + "umap_sample_37c_no_t.svg",
#     conda:
#         "envs/DB_AKC_scarf.yaml"
#     params:
#         threads=config["n_threads"],
#     script:
#         "src/smk/scarf_plot_37c_targets.py"


# rule seurat_explore_pca_loadings:
#     input:
#         seurat_object=rules.seurat_processing_all.output.seurat_object,
#     output:
#         csvs=...,
#         plots=...,
#     conda:
#         "envs/DB_AKC_scarf.yaml"
#     script:
#         "src/smk/seurat_explore_pca_loadings.R"


