import scarf
import pandas as pd


ds = scarf.DataStore(
    snakemake.input["zarr_ice"],
    default_assay="RNA",
    nthreads=snakemake.params["threads"],
)

ds.mark_hvgs(from_assay="RNA", top_n=1000)
ds.make_graph(feat_key="hvgs", return_ann_object=True)
ds.run_umap(from_assay="RNA")
ds.run_leiden_clustering(from_assay="RNA", resolution=.6)
