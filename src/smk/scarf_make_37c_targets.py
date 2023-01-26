import scarf
import pandas as pd

ds = scarf.DataStore(
    snakemake.input["zarr_37c"],
    default_assay="RNA",
    nthreads=snakemake.params["threads"],
)

############################################################
#  Prepare split of 37c samples between buffer_treatment:  #
############################################################

htos = ds.cells.fetch_all("hto")

hto_number = list(
    set(snakemake.config["sample_groups"]["37c"]).intersection(
        set(snakemake.config["sample_groups"]["with_triptolide"])
    )
)
if len(hto_number) > 1:
    raise Exception()
else:
    hto_number = str(int(hto_number[0]))
hto_37c_t = "HTO" + hto_number
is_37c_t = [True if x in hto_37c_t else False for x in htos]
ds.cells.insert("is_37c_t", is_37c_t, overwrite=True)

hto_number = list(
    set(snakemake.config["sample_groups"]["37c"]).intersection(
        set(snakemake.config["sample_groups"]["no_triptolide"])
    )
)
if len(hto_number) > 1:
    raise Exception()
else:
    hto_number = str(int(hto_number[0]))
hto_37c_no_t = "HTO" + hto_number
is_37c_no_t = [True if x in hto_37c_no_t else False for x in htos]
ds.cells.insert("is_37c_no_t", is_37c_no_t, overwrite=True)

#########################
#  Export 37c samples:  #
#########################

scarf.writers.SubsetZarr(
    in_zarr=snakemake.input["zarr_37c"],
    out_zarr=snakemake.output["zarr_37c_t"],
    cell_key="is_37c_t",
    reset_cell_filter=False,
    overwrite_existing_file=True,
).dump()

scarf.writers.SubsetZarr(
    in_zarr=snakemake.input["zarr_37c"],
    out_zarr=snakemake.output["zarr_37c_no_t"],
    cell_key="is_37c_no_t",
    reset_cell_filter=False,
    overwrite_existing_file=True,
).dump()

###########################
#  Process both samples:  #
###########################

ds = scarf.DataStore(
    snakemake.output["zarr_37c_t"],
    default_assay="RNA",
    nthreads=snakemake.params["threads"],
)

ds.mark_hvgs(from_assay="RNA", top_n=1000)
ds.make_graph(feat_key="hvgs", return_ann_object=True)
ds.run_umap(from_assay="RNA")
ds.run_leiden_clustering(from_assay="RNA", resolution=0.6)

ds = scarf.DataStore(
    snakemake.output["zarr_37c_no_t"],
    default_assay="RNA",
    nthreads=snakemake.params["threads"],
)

ds.mark_hvgs(from_assay="RNA", top_n=1000)
ds.make_graph(feat_key="hvgs", return_ann_object=True)
ds.run_umap(from_assay="RNA")
ds.run_leiden_clustering(from_assay="RNA", resolution=0.6)
