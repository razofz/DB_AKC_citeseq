# ---
# jupyter:
#   jupytext:
#     formats: ipynb,Rmd,R
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
# title: "Revision plots, IER manuscript"
# format: 
#   html:
#     code-fold: false
#     theme: journal
#     embed-resources: true
# jupyter: R
# author: "Rasmus Olofzon"
# toc: true
# toc-expand: true
# fig-cap-location: top
# ---

# ---

#
# ## Methods
#
# ### Differentially expressed genes (DEG)
#
# The samples were first compared cell cycle phase-wise, since the dataset showed clear separation based on phase classification, but its structure otherwise made biological sense. In the UMAP the two ice conditions co-located very clearly, which prompted the notion of treating them as one for the DEG analysis. This was corroborated by a DEG testing between the two ice conditions, which for the samples as wholes yielded three DEGs, non with significant p-values. Phase-wise no DEGs were found.
#
# Based off of that, three comparisons were carried out:
#
# 1. 37c _with_ triptolide VS ice
# 2. 37c _without_ triptolide VS ice
# 3. 37c _without_ triptolide VS 37c _without_ triptolide
#
# The comparisons were first made cell cycle phase-wise (G1, G2M, S), then the intersection of the found DEGs was taken. This would give the differentially expressed genes that most characterizes a given condition, without cell cycle-specific genes.
#
# The DEG testing was performed with Seurat's `FindMarkers` function, testing performed only on highly variable genes (HVGs). The intersection was taken with R's `intersect` function for sets.
#
# ### Gene signatures
#
# The gene signature scores were calculated and visualized similar to the procedure for the IER signature.
#
# ---

# # Setup

# + vscode={"languageId": "r"}
library(Seurat)
library(stringr)
library(dplyr)
library(ggplot2)

# + vscode={"languageId": "r"}
library(tidyr)
library(dplyr)
library(viridis)

# + vscode={"languageId": "r"}
library(patchwork)

# + vscode={"languageId": "r"}
sobj <- readRDS("../data/processed/seurat_object_w_stress_sig.rds")

# + vscode={"languageId": "r"}
# cols_features <- c("moccasin", "darkslategray")
cols_features <- c("lightgray", "red3")


# + vscode={"languageId": "r"}
phases <- list("G1", "G2M", "S")
names(phases) <- phases
phases
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
df <- table(sobj[[]][,c("sample", "Phase")]) %>% as.data.frame.matrix #%>% rbind(
# table(sobj[[]][,c("sample")]) %>% as.list
# )
df

# + vscode={"languageId": "r"}
df[["sample"]] <- table(sobj[["sample"]]) %>% as.list %>% as.numeric
df

# + vscode={"languageId": "r"}
df <- df %>% mutate(
    G1_frac = round(G1 / sample, digits = 2),
    G2M_frac = round(G2M / sample, digits = 2),
    S_frac = round(S / sample, digits = 2),
)
df

# + vscode={"languageId": "r"}
#| label: tbl-phase-sample
#| tbl-cap: "Cell cycle phase classification distribution per condition"
df <- table(sobj[[]][,c("sample", "Phase")]) %>% as.data.frame.matrix #%>% rbind(
# table(sobj[[]][,c("sample")]) %>% as.list
# )
df
df <- df %>% tibble::rownames_to_column(var = "sample") %>% pivot_longer(cols = c("G1", "G2M", "S"), names_to = "phase", values_to="n_cells")
df

# + vscode={"languageId": "r"}
ggplot(df, aes(fill=phase, y=n_cells, x=sample)) + 
    geom_bar(position="fill", stat="identity") +
    scale_fill_viridis(discrete = T, option = "cividis") +
    ylab("") +
    xlab("Condition")
ggsave(last_plot(), filename = "plots/cc_phase_distribution_per_condition.svg", device = "svg", units = "in", width = 6, height = 6)
# -

# # Compare ice conditions

# + vscode={"languageId": "r"}
WhichCells(sobj, expr = sample == "ice_t" & Phase == "S") %>% length

# + vscode={"languageId": "r"}
Idents(sobj) <- "sample"
DefaultAssay(sobj) <- "RNA"
degs <- FindMarkers(
    sobj,
    ident.1 = "ice_t",
    ident.2 = "ice_no_t",
    assay = "RNA",
    features = VariableFeatures(sobj)
) %>% arrange(desc(avg_log2FC))
degs
# -

# So, only three DEGs, with insignificant p-values. Should probably be a good argument for treating them as one for this.

# + vscode={"languageId": "r"}
get_degs_cc_wise <- \(idents = "sample", group_1 = "ice_t", group_2 = c("ice_t", "ice_no_t"), cc_phase = "S") {
    Idents(sobj) <- idents
    DefaultAssay(sobj) <- "RNA"

    # extract cell groups:
    cells_1 = WhichCells(sobj, expr = sample %in% group_1 & Phase == cc_phase)
    cells_2 = WhichCells(sobj, expr = sample %in% group_2 & Phase == cc_phase)
    foo_1 <- table(sobj[[]][cells_1,c("sample", "Phase")]) %>% as.data.frame.matrix
    foo_2 <- table(sobj[[]][cells_2,c("sample", "Phase")]) %>% as.data.frame.matrix

    # double-check that the correct cell groups are extracted,
    # will throw an error and fail if they are not correct:
    stopifnot(phase_distributions[group_2, cc_phase] == foo_2[group_2, cc_phase])

    degs <- FindMarkers(
        sobj[["RNA"]],
        cells.1 = cells_1,
        cells.2 = cells_2,
        features = VariableFeatures(sobj)
    ) %>% arrange(desc(avg_log2FC))
    degs %>% return
}

# + vscode={"languageId": "r"}
get_degs_cc_wise(group_1 = "ice_t", group_2 = "ice_no_t", cc_phase = "S")
get_degs_cc_wise(group_1 = "ice_t", group_2 = "ice_no_t", cc_phase = "G1")
get_degs_cc_wise(group_1 = "ice_t", group_2 = "ice_no_t", cc_phase = "G2M")
# -

# So yeah, no DEGs at all doing it phase-wise.

# # Compare 37c conditions

# + vscode={"languageId": "r"}
sobj[[]][["sample"]] %>% unique
# -

# - 37c_t vs ice
# - 37c_no_t vs ice
# - 37c_t vs 37c_no_t

# + vscode={"languageId": "r"}
ice <- c("ice_t", "ice_no_t")
comparisons <- list(
    "37c_t" = list("37c_t", ice),
    "37c_no_t" = list("37c_no_t", ice),
    "37c_t_vs_no_t" = list("37c_t", "37c_no_t")
)
# degs <- list(
#     "37c_t" = "",
#     "37c_no_t" = "",
#     "37c_t_vs_no_t" = ""
# )
comparisons
# degs

# + vscode={"languageId": "r"}
get_degs_cc_wise(group_1 = comparisons[["37c_t"]][[1]], group_2 = comparisons[["37c_t"]][[2]], cc_phase = "S") %>% head

# + vscode={"languageId": "r"}
phase_distributions <- table(sobj[[]][,c("sample", "Phase")]) %>% as.data.frame.matrix #%>% rbind(
phase_distributions

# + vscode={"languageId": "r"}
degs <- lapply(comparisons, FUN = \(comp) {
    lapply(phases, FUN = \(phase) {
        group_1 <- comp[[1]]
        group_2 <- comp[[2]]
        print(group_1)
        print(group_2)
        print(phase)
        cells_1 = WhichCells(sobj, expr = sample %in% group_1 & Phase == phase)
        cells_2 = WhichCells(sobj, expr = sample %in% group_2 & Phase == phase)
        print(length(cells_1))
        print(length(cells_2))
        foo_1 <- table(sobj[[]][cells_1,c("sample", "Phase")]) %>% as.data.frame.matrix
        foo_2 <- table(sobj[[]][cells_2,c("sample", "Phase")]) %>% as.data.frame.matrix
        print(foo_1)
        print(foo_2)
        stopifnot(phase_distributions[group_2, phase] == foo_2[group_2, phase])
        print("------")
    })
})

# + vscode={"languageId": "r"}
degs <- lapply(comparisons, FUN = \(comp) {
    lapply(phases, FUN = \(phase) {
        get_degs_cc_wise(group_1 = comp[[1]], group_2 = comp[[2]], cc_phase = phase)
    })
})

# + vscode={"languageId": "r"}
names(degs)
names(degs[["37c_t"]])

# + vscode={"languageId": "r"}
degs[["37c_t"]][["S"]] %>% head
degs[["37c_t_vs_no_t"]][["S"]] %>% head
# -

# ## Save DEG tables

# + vscode={"languageId": "r"}
write.csv(degs[["37c_t"]][["S"]], file = "plots/1.tsv")

# + vscode={"languageId": "r"}
lapply(names(comparisons), FUN = \(comp) {
    comp_name <- comp
    comp <- comparisons[[comp_name]]
    if (length(comp[[2]] > 1)) comp[[2]] <- str_flatten(comp[[2]], "_")
    # print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>")
    # print(comp_name)
    # print(comp)
    lapply(phases, FUN = \(phase) {
        # print(comp)
        # print(phase)
        filename <- str_c(
            "plots/", comp[[1]], "_VS_", comp[[2]], "_phase_", phase, ".tsv"
        )
        print(filename)
        # print(dim(degs[[comp_name]][[phase]]))
        # print(head(degs[[comp_name]][[phase]]))
        # print("------")
        write.csv(degs[[comp_name]][[phase]], file = filename)
        # return(file.exists(filename))
        return(dim(degs[[comp_name]][[phase]]))
    })
})
# -

# ## Overlapping DEGs (intersection between cell cycle phases)

# + vscode={"languageId": "r"}
comparisons

# + vscode={"languageId": "r"}
tmp <- lapply(degs[["37c_t"]], rownames)
lapply(tmp, length)
intersect(tmp[["G1"]], tmp[["G2M"]]) %>% length
intersect(tmp[["G1"]], tmp[["S"]]) %>% length
intersect(tmp[["S"]], tmp[["G2M"]]) %>% length

a <- intersect(tmp[["G1"]], tmp[["G2M"]])
b <- intersect(tmp[["G1"]], tmp[["S"]])
c <- intersect(tmp[["S"]], tmp[["G2M"]])

intersect(a, intersect(b, c)) %>% length
intersect(a, intersect(b, c))

# + vscode={"languageId": "r"}
# degs[["37c_t"]][["S"]][ intersect(a, intersect(b, c)), degs[["37c_t"]][["S"]]$avg_log2FC > 0]
# degs[["37c_t"]][["S"]][ intersect(a, intersect(b, c)) & "avg_log2FC" < 0, ]
degs[["37c_t"]][["S"]][ intersect(a, intersect(b, c)), ] %>% dim
degs[["37c_t"]][["S"]][ intersect(a, intersect(b, c)), ] %>% head

# + vscode={"languageId": "r"}
overlapping_degs <- \(comp = "37c_t") {
    tmp <- lapply(degs[[comp]], rownames)
    overlaps <- intersect( tmp[["G1"]], intersect(tmp[["G2M"]], tmp[["S"]]))
    return(degs[[comp]][["G1"]])
}

# + vscode={"languageId": "r"}
overlaps <- lapply(names(comparisons), FUN = \(comp) overlapping_degs(comp = comp))
names(overlaps) <- names(comparisons)
lapply(overlaps, dim)
lapply(overlaps, head)

# + vscode={"languageId": "r"}
lapply(names(comparisons), FUN = \(comp) {
    write.csv(overlaps[[comp]], file = str_c(
        "plots/", "deg_overlap_", comp, ".tsv"
    ))
})
# -

# ---

# + vscode={"languageId": "r"}
feats <- c(
    degs[["37c_t"]][["G1"]] %>% head(., n = 15) %>% rownames,
    degs[["37c_t"]][["G1"]] %>% tail(., n = 15) %>% rownames
)
p1 <- DotPlot(sobj, features = feats, cols = cols_features) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
feats <- c(
    degs[["37c_t"]][["G2M"]] %>% head(., n = 15) %>% rownames,
    degs[["37c_t"]][["G2M"]] %>% tail(., n = 15) %>% rownames
)
p2 <- DotPlot(sobj, features = feats, cols = cols_features) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
feats <- c(
    degs[["37c_t"]][["S"]] %>% head(., n = 15) %>% rownames,
    degs[["37c_t"]][["S"]] %>% tail(., n = 15) %>% rownames
)
p3 <- DotPlot(sobj, features = feats, cols = cols_features) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# + vscode={"languageId": "r"}
plot <- (p1 + grid::textGrob("G1")) / p2 / p3
plot

# + vscode={"languageId": "r"}
ggsave(plot, filename = "plots/1.svg", device = "svg", units = "in", width = 12, height = 12)

# + vscode={"languageId": "r"}
one_dot_plot <- \(comparison = "37c_t", n_genes = 15, phase = "S") {
    feats <- c(
        degs[[comparison]][[phase]] %>% head(., n = n_genes) %>% rownames,
        degs[[comparison]][[phase]] %>% tail(., n = n_genes) %>% rownames
    )
    plot <- DotPlot(sobj, features = feats, cols = cols_features) +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
        ) +
        labs(
            # title = str_c("Comparison ", comparison, ":"),
            subtitle = str_c("Cell cycle phase ", phase, ":"),
            y = "",
            x = ""
            # y = "Condition",
            # x = "DEGs"
        )
    return(plot)
}
one_dot_plot()

# + vscode={"languageId": "r"}
comp <- "37c_t"
n_genes <- 15
plot <- one_dot_plot(comparison = comp, phase = "G1") /
one_dot_plot(comparison = comp, phase = "G2M") /
one_dot_plot(comparison = comp, phase = "S") + plot_annotation(
  title = str_c(
    "DEGs between ",
    comparisons[[comp]][[1]],
    " and ",
    "(", str_flatten(comparisons[[comp]][[2]], " + "), ")"
  ),
  subtitle = "",
  caption = str_c(
    "The top ", n_genes, " DEGs for ", comp[[1]],
    " and the top ", n_genes, " DEGs for ", 
    str_flatten(comparisons[[comp]][[2]], " + ")
  )
)
plot

# + vscode={"languageId": "r"}
ggsave(plot, filename = "plots/1.svg", device = "svg", units = "in", width = 12, height = 12)

# + vscode={"languageId": "r"}
all_phases_dot_plot <- \(comparison = "37c_t", n_genes = 15) {
  plot <- one_dot_plot(comparison = comparison, phase = "G1") /
  one_dot_plot(comparison = comparison, phase = "G2M") /
  one_dot_plot(comparison = comparison, phase = "S") + plot_annotation(
    title = str_c(
      "DEGs between ",
      comparisons[[comparison]][[1]],
      " and ",
      str_flatten(comparisons[[comparison]][[2]], " + ")
    ),
    subtitle = "Cell cycle phase-wise",
    caption = str_c(
      "The top ", n_genes, " DEGs for ", comparisons[[comparison]][[1]],
      " (left) and the top ", n_genes, " DEGs for ", 
      str_flatten(comparisons[[comparison]][[2]], " + "),
      " (right)"
    )
  )
  plot %>% return
}
all_phases_dot_plot()
# -

# ---
#
# # Collated dot plots

# + vscode={"languageId": "r"}
plot <- all_phases_dot_plot(comparison = "37c_t")
plot
ggsave(plot, filename = "plots/37c_t.svg", device = "svg", units = "in", width = 12, height = 12)

# + vscode={"languageId": "r"}
plot <- all_phases_dot_plot(comparison = "37c_no_t")
plot
ggsave(plot, filename = "plots/37c_no_t.svg", device = "svg", units = "in", width = 12, height = 12)

# + vscode={"languageId": "r"}
plot <- all_phases_dot_plot(comparison = "37c_t_vs_no_t")
plot
ggsave(plot, filename = "plots/37c_t_vs_no_t.svg", device = "svg", units = "in", width = 12, height = 12)
# -

# ## Dot plots overlapping DEGs

# + vscode={"languageId": "r"}
one_dot_plot_overlaps <- \(comparison = "37c_t", n_genes = 15) {
    feats <- c(
        overlaps[[comparison]] %>% head(., n = n_genes) %>% rownames,
        overlaps[[comparison]] %>% tail(., n = n_genes) %>% rownames
    )
    plot <- DotPlot(sobj, features = feats, cols = cols_features) +
        theme(
            axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)
        ) +
        labs(
            title = str_c(
                "DEGs between ",
                comparisons[[comparison]][[1]],
                " and ",
                str_flatten(comparisons[[comparison]][[2]], " + ")
            ),
            subtitle = str_c(
                "Overlaps/intersection between cell cycle phase-wise DEG analyses"
            ),
            caption = str_c(
                "The top ", n_genes, " DEGs for ", comparisons[[comparison]][[1]],
                " (left) and the top ", n_genes, " DEGs for ", 
                str_flatten(comparisons[[comparison]][[2]], " + "),
                " (right)"
            ),
            y = "Condition",
            x = "DEGs"
        )
    return(plot)
}
one_dot_plot_overlaps()

# + vscode={"languageId": "r"}
names(comparisons)

# + vscode={"languageId": "r"}
one_dot_plot_overlaps("37c_t")
one_dot_plot_overlaps("37c_t_vs_no_t")
one_dot_plot_overlaps("37c_no_t")

# + vscode={"languageId": "r"}
lapply(names(overlaps), FUN = \(comp) {
    ggsave(
        one_dot_plot_overlaps(comparison = comp),
        filename = str_c("plots/", "deg_overlap_", comp, ".svg"),
        device = "svg", units = "in", width = 12, height = 5
    )
})
# -

# ---
#
# # QC plots

# + vscode={"languageId": "r"}
sobj[[]] %>% colnames

# + vscode={"languageId": "r"}
feats_of_interest <- list(
  "nCount_RNA", "nFeature_RNA", "nCount_ADT", "nCount_HTO", "percent.mt", "stress_signature1"
)
feats_of_interest

# + vscode={"languageId": "r"}
qc_plot <- \(feat = "nCount_RNA") {
    VlnPlot(sobj, features = feat) +
        scale_fill_viridis(discrete = T, option = "plasma")
}
qc_plot()

# + vscode={"languageId": "r"}
plots <- lapply(feats_of_interest, qc_plot)

# + vscode={"languageId": "r"}
names(plots) <- feats_of_interest

# + vscode={"languageId": "r"}
lapply(feats_of_interest, FUN = \(feat) ggsave(plots[[feat]],
    filename = str_c("plots/qc_", feat, ".svg"),
    device = "svg",
    units = "in",
    width = 6,
    height = 6)
)

# + vscode={"languageId": "r"}
plots
# -

# ---
#
# # Signatures

# + vscode={"languageId": "r"}
signature_file <- "../data/raw/Genesets.gmx.txt"

# + vscode={"languageId": "r"}
dirname(signature_file) %>% dir

# + vscode={"languageId": "r"}
sigs_dirty <- read.delim(signature_file)
sigs_dirty %>% head

# + vscode={"languageId": "r"}
sigs <- list()

# + vscode={"languageId": "r"}
sigs[["cc"]] <- sigs_dirty[[colnames(sigs_dirty)[
    str_detect("cycle", string = colnames(sigs_dirty))]
]]
sigs[["young"]] <- sigs_dirty[[colnames(sigs_dirty)[
    str_detect("Young", string = colnames(sigs_dirty))]
]]
sigs[["aged"]] <- sigs_dirty[[colnames(sigs_dirty)[
    str_detect("Aged", string = colnames(sigs_dirty))]
]]
sigs[["stress"]] <- sigs_dirty[[colnames(sigs_dirty)[
    str_detect("Stress", string = colnames(sigs_dirty))]
]]
sigs

# + vscode={"languageId": "r"}
saved_names <- names(sigs)
sigs <- lapply(sigs, FUN = \(s) s[s %>% str_length > 0 & s != "na"])
# sigs <- lapply(names(sigs), FUN = \(s) sigs[[s]][sigs[[s]] %>% str_length > 0 & sigs[[s]] != "na"])
# names(sigs) <- saved_names
sigs

# + vscode={"languageId": "r"}
lapply(sigs, FUN = \(s) s %in% rownames(sobj) %>% table)
# lapply(sigs, FUN = \(s) lapply(s,
#     FUN = \(gene) str_detect(pattern = gene, string = rownames(sobj)) %>% sum
#     ) %>% table
# )

# + vscode={"languageId": "r"}
s <- "cc"
new <- AddModuleScore(sobj, features = sigs, name = names(sigs))
new

# + vscode={"languageId": "r"}
new[[]] %>% head

# + vscode={"languageId": "r"}
cols <- c("gray", "lightcoral", "red")


# + vscode={"languageId": "r"}
cc <- FeaturePlot(new, features = "cc1", coord.fixed = T, order = T, cols = cols) +
    ggtitle("Cell cycle signature")
young <- FeaturePlot(new, features = "young2", coord.fixed = T, order = T, cols = cols) +
    ggtitle("Young signature")
aged <- FeaturePlot(new, features = "aged3", coord.fixed = T, order = T, cols = cols) +
    ggtitle("Aged signature")
cc
young
aged

# + vscode={"languageId": "r"}
ggsave(cc, filename = str_c("plots/signature_", "cc", ".svg"),
    device = "svg", units = "in", width = 6, height = 6)
ggsave(young, filename = str_c("plots/signature_", "young", ".svg"),
    device = "svg", units = "in", width = 6, height = 6)
ggsave(aged, filename = str_c("plots/signature_", "aged", ".svg"),
    device = "svg", units = "in", width = 6, height = 6)
