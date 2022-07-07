import scarf
import pandas as pd


ds = scarf.DataStore(
    snakemake.input["zarr_ice"],
    default_assay="RNA",
    nthreads=snakemake.params["threads"],
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='RNA_leiden_cluster',
    savename=snakemake.output['plot_umap_clusters'],
    show_fig=False
)

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='Phase',
    savename=snakemake.output['plot_umap_phase'],
    show_fig=False
)

# ds.plot_layout(
#     layout_key='RNA_UMAP',
#     color_by='buffer_treatment',
#     savename=snakemake.output['plot_umap_buffer_treatment'],
#     show_fig=False
# )

# ds.plot_layout(
#     layout_key='RNA_UMAP',
#     color_by='hto',
#     savename=snakemake.output['plot_umap_hto'],
#     show_fig=False
# )

ds.plot_layout(
    layout_key='RNA_UMAP',
    color_by='sample',
    savename=snakemake.output['plot_umap_sample'],
    show_fig=False
)
