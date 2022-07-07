import scarf
import pandas as pd

metadata = pd.read_csv(
    snakemake.input["seurat_metadata"], sep="\t", index_col=0
)

reader = scarf.CrH5Reader(snakemake.input["filtered_h5"])
# reader = scarf.CrH5Reader(
#     "../../raw/cellranger-mm10/filtered_feature_bc_matrix.h5"
# )

writer = scarf.CrToZarr(
    reader,
    zarr_fn=snakemake.output["zarr_all"],
    # zarr_fn="ice.zarr",
    chunk_size=(2000, 1000),
)
writer.dump(batch_size=1000)

ds = scarf.DataStore(
    snakemake.output["zarr_all"],
    # "ice.zarr",
    default_assay="RNA",
    nthreads=snakemake.params["threads"],
)

# do same cell filtering here in Scarf as was done in Seurat:
md = ds.cells.to_pandas_dataframe(columns=ds.cells.columns)
cell_in_seurat_metadata = md.names.isin(metadata.index)
ds.cells.update_key(cell_in_seurat_metadata, "I")

# now we can insert metadata from Seurat directly:
ds.cells.insert("hto", values=metadata.hto)
ds.cells.insert("sample", values=metadata["sample"])
ds.cells.insert("sample_longname", values=metadata.sample_longname)
ds.cells.insert("buffer_treatment", values=metadata.buffer_treatment)
ds.cells.insert("incubation_method", values=metadata.incubation_method)
ds.cells.insert("Phase", values=metadata.Phase)
ds.cells.insert("S.Score", values=metadata["S.Score"])
ds.cells.insert("G2M.Score", values=metadata["G2M.Score"])
ds.cells.insert(
    "seurat_clusters", values=metadata.seurat_clusters, fill_value=-1
)

# prepare split of dataset between incubation methods:
htos = ds.cells.fetch_all("hto")

htos_in_37c = ["HTO" + str(i) for i in snakemake.config["sample_groups"]["37c"]]
is_37c = [True if x in htos_in_37c else False for x in htos]
# is_37c = [True if x in ["HTO2", "HTO4"] else False for x in htos]
ds.cells.insert("is_37c", is_37c, overwrite=True)

htos_in_ice = ["HTO" + str(i) for i in snakemake.config["sample_groups"]["ice"]]
is_ice = [True if x in htos_in_ice else False for x in htos]
# is_ice = [True if x in ["HTO1", "HTO3"] else False for x in htos]
ds.cells.insert("is_ice", is_ice, key="I", overwrite=True)

# export ice cells:
scarf.writers.SubsetZarr(
    in_zarr=snakemake.output["zarr_all"],
    out_zarr=snakemake.output["zarr_ice"],
    cell_key="is_ice",
    reset_cell_filter=False,
    overwrite_existing_file=True,
)

# export 37c cells:
scarf.writers.SubsetZarr(
    in_zarr=snakemake.output["zarr_all"],
    out_zarr=snakemake.output["zarr_37c"],
    cell_key="is_37c",
    reset_cell_filter=False,
    overwrite_existing_file=True,
)

# ds.mark_hvgs()
# ds.make_graph(feat_key="hvgs")
# ds.run_umap()
# ds.run_leiden_clustering(resolution=1)
