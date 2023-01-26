import scarf
import pandas as pd


########################
#  Plot 37c_t sample:  #
########################

ds = scarf.DataStore(
    snakemake.input["zarr_37c_t"],
    default_assay="RNA",
    nthreads=snakemake.params["threads"],
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='RNA_leiden_cluster',
    savename=snakemake.output['plot_umap_clusters_37c_t'],
    show_fig=False
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='Phase',
    savename=snakemake.output['plot_umap_phase_37c_t'],
    show_fig=False
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='sample',
    savename=snakemake.output['plot_umap_sample_37c_t'],
    show_fig=False
)

###########################
#  Plot 37c_no_t sample:  #
###########################

ds = scarf.DataStore(
    snakemake.input["zarr_37c_no_t"],
    default_assay="RNA",
    nthreads=snakemake.params["threads"],
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='RNA_leiden_cluster',
    savename=snakemake.output['plot_umap_clusters_37c_no_t'],
    show_fig=False
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='Phase',
    savename=snakemake.output['plot_umap_phase_37c_no_t'],
    show_fig=False
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='sample',
    savename=snakemake.output['plot_umap_sample_37c_no_t'],
    show_fig=False
)
