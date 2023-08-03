# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# ---
# title: "IER DEGs"
# format: 
#   html:
#     code-fold: true
#     theme: journal
#     embed-resources: true
#   pdf:
#     papersize: letter
# jupyter: R
# author: "Rasmus Olofzon"
# toc: true
# toc-expand: true
# fig-cap-location: top
# ---

# # Summary and overview

# Highlights:
#
# - UMAPs with sample and cell cycle phase visualized: @fig-dimplots
# - Metadata for dataset: @tbl-metadata
# - Cell cycle phase classification distribution per condition: @tbl-phase-sample
# - Previously performed DEG analysis, comparisons shown here explained in this diagram: @fig-diagram
#     - For **G1** phase: @fig-deg-old-g1
#     - For **G2M** phase: @fig-deg-old-g2m
#     - For **S** phase: @fig-deg-old-s
# - New DEG analysis (explanatory diagram in @fig-diagram-new)
#     - Table with DEGs: @tbl-new-degs
#     - Figure with top DEGs shown: @fig-deg-new

# # Setup

# + vscode={"languageId": "r"}
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)

# + vscode={"languageId": "r"}
sobj <- readRDS("../data/processed/seurat_object_w_stress_sig.rds")

# + vscode={"languageId": "r"}
sobj

# + vscode={"languageId": "r"}
#| echo: false
# cond <- list("37c_t", "ice_t")
# deg_files <- lapply(list("37c_v_ice" = cond, "ice_v_37c" = rev(cond)), FUN = \(cond) {
#     phases <- list("G1", "G2M", "S")
#     names(phases) <- phases
#     lapply(phases, FUN = \(x) str_c("../data/processed/deg/deg_", x, "_", cond[[1]], "_vs_", cond[[2]], ".tsv"))
# })
# deg_files

# + vscode={"languageId": "r"}
#| echo: false
# degs <- lapply(list("37c_v_ice" = "37c_v_ice", "ice_vs_37c" = "ice_vs_37c"), FUN = \(comparison) {
#     phases <- list("G1", "G2M", "S")
#     names(phases) <- phases
#     lapply(phases, FUN = \(phase) read.table(deg_files[[comparison]][[phase]]))
# })

# + vscode={"languageId": "r"}
cols_features <- c("moccasin", "darkslategray")

# + vscode={"languageId": "r"}
phases <- list("G1", "G2M", "S")
names(phases) <- phases
deg_files <- lapply(phases, FUN = \(x) str_c("../data/processed/deg/deg_", x, "_ice_t_vs_37c_t.tsv"))
# deg_files

# + vscode={"languageId": "r"}
degs <- lapply(phases, FUN = \(phase) read.table(deg_files[[phase]]))
# -

# ---
#
# # Visualization of dataset with conditions, for reference
#
# Metadata:

# + vscode={"languageId": "r"}
#| label: tbl-metadata
#| tbl-cap: "Metadata for dataset"
sobj[[]] %>% head

# + vscode={"languageId": "r"}
#| label: fig-dimplots
#| fig-cap: "UMAP coloured on metadata"
#| fig-subcap: 
#|   - "Coloured on condition/sample"
#|   - "Coloured on cell cycle phases"
#| layout-ncol: 2
DimPlot(sobj, group.by = "sample") + coord_fixed()
DimPlot(sobj, group.by = "Phase") + coord_fixed()

# + vscode={"languageId": "r"}
#| label: tbl-phase-sample
#| tbl-cap: "Cell cycle phase classification distribution per condition"
table(sobj[[]][,c("sample", "Phase")]) %>% as.data.frame.matrix
# -

# ---
#
# # Old DEG analysis
#
# From former analysis we performed differential gene expression testing cell cycle phase-wise. The comparison was between conditions/samples. This yielded the following set of files:
#

# + vscode={"languageId": "r"}
system("ls ../data/processed/deg/", intern = T)
# -

# ```{mermaid}
# %%| label: fig-diagram
# %%| fig-cap: "Schema for (old) DEG comparisons to follow."
# flowchart TB
#   subgraph ice with triptolide
#   G1_ice[G1]
#   G2M_ice[G2M]
#   S_ice[S]
#   end
#   subgraph 37c with triptolide
#   G1_37c[G1] --- G1_ice
#   G2M_37c[G2M] --- G2M_ice
#   S_37c[S] --- S_ice
#   end
# ```

# For each file, both up- and down-regulated genes are included. Here we have loaded the comparison between 37c with triptolide and ice with triptolide (see @fig-diagram). So first we have the top 6 ice-specific, for the G1 phase:

# + vscode={"languageId": "r"}
#| code-fold: false
degs[["G1"]] %>% head
# -

# And the top 6 37c-specific, for the G1 phase as well:

# + vscode={"languageId": "r"}
#| code-fold: false
degs[["G1"]] %>% tail
# -

# ---
#
# ## DEGs visualized on UMAP
#
# ### G1

# + vscode={"languageId": "r"}
#| label: fig-deg-old-g1
#| fig-cap: "Top DEGs for phase **G1**."
#| fig-subcap: 
#|   - "Top 6 DEGs specific to **ice w/ t** condition"
#|   - "Top 6 DEGs specific to **37c w/ t** condition"
#| layout-ncol: 1
FeaturePlot(sobj, features = degs[["G1"]] %>% head %>% rownames, coord.fixed = T, order = T, cols = cols_features)
FeaturePlot(sobj, features = degs[["G1"]] %>% tail %>% rownames, coord.fixed = T, order = T, cols = cols_features)
# -

# ### G2M

# + vscode={"languageId": "r"}
#| label: fig-deg-old-g2m
#| fig-cap: "Top DEGs for phase **G2M**."
#| fig-subcap: 
#|   - "Top 6 DEGs specific to **ice w/ t** condition"
#|   - "Top 6 DEGs specific to **37c w/ t** condition"
#| layout-ncol: 1
FeaturePlot(sobj, features = degs[["G2M"]] %>% head %>% rownames, coord.fixed = T, order = T, cols = cols_features)
FeaturePlot(sobj, features = degs[["G2M"]] %>% tail %>% rownames, coord.fixed = T, order = T, cols = cols_features)
# -

# ### S

# + vscode={"languageId": "r"}
#| label: fig-deg-old-s
#| fig-cap: "Top DEGs for phase **S**."
#| fig-subcap: 
#|   - "Top 6 DEGs specific to **ice w/ t** condition"
#|   - "Top 6 DEGs specific to **37c w/ t** condition"
#| layout-ncol: 1
FeaturePlot(sobj, features = degs[["S"]] %>% head %>% rownames, coord.fixed = T, order = T, cols = cols_features)
FeaturePlot(sobj, features = degs[["S"]] %>% tail %>% rownames, coord.fixed = T, order = T, cols = cols_features)
# -

# ---
#
# # New DEG analysis
#
# Specifically done _not_ cell cycle phase-wise. Will be done for:
#
# - 37c with tript _VS_ ice with tript

# ```{mermaid}
# %%| label: fig-diagram-new
# %%| fig-cap: "Schema for new DEG comparison."
# flowchart TB
#   subgraph ice with triptolide
#   ice[G1, G2M, S]
#   end
#   subgraph 37c with triptolide
#   37c[G1, G2M, S] --- ice
#   end
# ```

# ##  37c with tript _vs_ ice with tript

# Sample names used in object:

# + vscode={"languageId": "r"}
sobj[[]][["sample"]] %>% unique

# + vscode={"languageId": "r"}
#| code-fold: false
Idents(sobj) <- "sample"
DefaultAssay(sobj) <- "RNA"
new_degs <- FindMarkers(sobj, ident.1 = "37c_t", ident.2 = "ice_t", assay = "RNA", features = VariableFeatures(sobj)) %>% arrange(desc(avg_log2FC))
# -

# ### DEG table

# + vscode={"languageId": "r"}
#| label: tbl-new-degs
#| tbl-cap: "Differentially expressed genes between `37c_t` and `ice_t`. Positive fold change = up in `37c_t`, negative fold change = up in `ice_t`."
new_degs
# -

# ### DEG visualizations

# + vscode={"languageId": "r"}
#| label: fig-deg-new
#| fig-cap: "Top DEGs between 37c with triptolide and ice with triptolide."
#| fig-subcap: 
#|   - "Top 6 DEGs specific to **`37c_t`** condition"
#|   - "Top 6 DEGs specific to **`ice_t`** condition"
#| layout-ncol: 1
FeaturePlot(sobj, features = new_degs %>% head %>% rownames, coord.fixed = T, order = T, cols = cols_features)
FeaturePlot(sobj, features = new_degs %>% tail %>% rownames, coord.fixed = T, order = T, cols = cols_features)
