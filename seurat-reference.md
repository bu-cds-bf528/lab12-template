# Seurat Cheatsheet

Functions and patterns from the Seurat practice lab and the basic workflow.
Replace `obj` with your object's name.

## Load data

Load each sample into its own object. The `project` name becomes the sample
label in the `orig.ident` metadata column.

```r
counts1 <- Read10X(data.dir = "path/to/sample1/filtered_feature_bc_matrix")
sample1 <- CreateSeuratObject(counts = counts1, project = "sample1")

counts2 <- Read10X(data.dir = "path/to/sample2/filtered_feature_bc_matrix")
sample2 <- CreateSeuratObject(counts = counts2, project = "sample2")
```

Run each QC step on each sample separately, for example:

```r
sample1[["percent.mt"]] <- PercentageFeatureSet(sample1, pattern = "^mt-")
sample2[["percent.mt"]] <- PercentageFeatureSet(sample2, pattern = "^mt-")
```

## Look inside the object

| Function | Returns |
| --- | --- |
| `obj` | A short summary |
| `ncol(obj)`, `nrow(obj)` | Number of cells, number of genes |
| `Cells(obj)`, `Features(obj)` | Cell barcodes, gene names |
| `Assays(obj)`, `Layers(obj)` | Assays; layers (`counts`, `data`, `scale.data`) |
| `obj@meta.data` | Metadata: one row per cell |
| `Reductions(obj)` | Stored reductions (PCA, UMAP) |
| `Idents(obj)` | Each cell's current identity |

## Metadata

```r
obj$stim                      # one column, as a vector
obj[["stim"]]                 # one column, as a data frame
obj[[col_name, drop = TRUE]]  # column name stored in a variable, as a vector
```

Use `$` for everyday reading and adding columns. Use `[[ ]]` when the column
name is in a variable.

```r
table(obj$orig.ident)                     # cells per sample
table(Idents(obj), obj$orig.ident)        # cells per cluster per sample
prop.table(table(Idents(obj), obj$flag), margin = 1)   # fractions within each cluster
```

### Add columns

```r
obj$new_column <- values                                    # one value per cell
obj$is_type    <- obj$annotation %in% c("Type A", "Type B")  # TRUE/FALSE
obj$group      <- paste(obj$annotation, obj$orig.ident, sep = "_")
```

### Percentage of counts from a gene family

```r
grep("^mt-", rownames(obj), value = TRUE)          # always check what matched
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^mt-")
```

- Gene names are case-sensitive: mouse mitochondrial genes are `mt-`, human are
  `MT-`. The wrong case matches nothing and gives 0% for every cell, with no
  error.
- A loose pattern can match unrelated genes (for example `^MT` matches
  `MTOR`). Check the matches before using them.
- You can also pass an explicit gene list: `PercentageFeatureSet(obj, features = my_genes)`.

## Subset

```r
x <- subset(obj, subset = orig.ident == "sample1")
x <- subset(obj, subset = annotation == "Type A" & nFeature_RNA > 1000)
x <- subset(obj, subset = annotation %in% c("Type A", "Type B"))
x <- subset(obj, idents = c("0", "3"))            # by current identity
ncol(x)                                           # confirm the cell count
```

### Flag cells instead of removing them

Flagging keeps every cell, so you can see what a cutoff *would* remove.

```r
obj$fails_cutoff <- obj$nFeature_RNA < 500        # TRUE = would be removed
table(obj$fails_cutoff)
table(Idents(obj), obj$fails_cutoff)
```

A per-cluster summary table:

```r
library(dplyr)
data.frame(cluster = Idents(obj), flagged = obj$fails_cutoff) |>
  group_by(cluster) |>
  summarise(total = n(), flagged = sum(flagged),
            fraction_flagged = round(flagged / total, 2)) |>
  arrange(desc(fraction_flagged))
```

## Merge samples

```r
obj <- merge(sample1, y = sample2, add.cell.ids = c("sample1", "sample2"))
obj <- JoinLayers(obj)   # combine per-sample layers into one
table(obj$orig.ident)    # confirm both samples are there
```

`add.cell.ids` keeps barcodes from different samples from colliding.

## The standard workflow

```r
obj <- NormalizeData(obj)
obj <- FindVariableFeatures(obj, nfeatures = 2000)
obj <- ScaleData(obj)                       # required before PCA
obj <- RunPCA(obj, npcs = 50)
ElbowPlot(obj, ndims = 50)

obj <- FindNeighbors(obj, dims = 1:30)
obj <- FindClusters(obj, resolution = c(0.2, 0.4, 0.6, 0.8, 1.0))
obj <- RunUMAP(obj, dims = 1:30)
```

### Choose a resolution

```r
library(clustree)
clustree(obj, prefix = "RNA_snn_res.")
Idents(obj) <- "RNA_snn_res.0.6"            # set your chosen resolution
```

`FindClusters()` with several resolutions sets the identities to the *last*
one. Set them yourself afterward.

## Plots

```r
# QC
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        group.by = "orig.ident", pt.size = 0)
FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# UMAPs
DimPlot(obj, label = TRUE)                           # current identities
DimPlot(obj, group.by = "orig.ident")                # color by a metadata column
DimPlot(obj, split.by = "orig.ident", label = TRUE)  # one panel per sample

# A gene or metadata column on the UMAP, with cluster numbers
FeaturePlot(obj, features = "nFeature_RNA", label = TRUE, order = TRUE)
FeaturePlot(obj, features = "fails_cutoff", label = TRUE, order = TRUE)

# Compare a metric between flagged and kept cells
VlnPlot(obj, features = c("nCount_RNA", "percent.mt"),
        group.by = "fails_cutoff", pt.size = 0)
```

If `FeaturePlot()` errors on a TRUE/FALSE column, convert it first:
`obj$flag_num <- as.numeric(obj$fails_cutoff)`.


# Gene expression, one column per group (e.g. flagged vs. kept)
`FeaturePlot(obj, features = c("GeneA", "GeneB"), split.by = "fails_cutoff", order = TRUE)`


### Draw cutoffs on a plot

Plot one metric at a time, so each line sits on the right panel:

```r
VlnPlot(obj, features = "nFeature_RNA", pt.size = 0) +
  geom_hline(yintercept = c(500, 5000), linetype = "dashed", colour = "red")

VlnPlot(obj, features = "percent.mt", pt.size = 0) +
  geom_hline(yintercept = 10, linetype = "dashed", colour = "red")
```

Add `group.by = "orig.ident"` to show one violin per sample instead of per
cluster.

### Custom plots with ggplot

`FetchData()` pulls metadata and gene expression into a data frame:

```r
library(ggplot2)
df <- FetchData(obj, vars = c("nFeature_RNA", "percent.mt", "orig.ident"))
df$cluster <- Idents(obj)

ggplot(df, aes(nFeature_RNA)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = c(500, 5000), linetype = "dashed", colour = "red") +
  facet_wrap(~ orig.ident) +
  theme_classic()
```

## Markers and annotation

```r
markers <- FindAllMarkers(obj, min.pct = 0.25, logfc.threshold = 0.5)

top_markers <- markers |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 5)

DotPlot(obj, features = c("GeneA", "GeneB", "GeneC")) + RotatedAxis()
FeaturePlot(obj, features = c("GeneA", "GeneB"))
```

Top markers for selected clusters only:

```r
my_clusters <- c("2", "5", "7")

markers |>
  filter(cluster %in% my_clusters) |>
  group_by(cluster) |>
  slice_max(avg_log2FC, n = 5)
```

### Name clusters

```r
new_names <- c("0" = "Label A", "1" = "Label B")   # one per cluster
obj <- RenameIdents(obj, new_names)
obj$my_labels <- Idents(obj)
```

### Optional: reference-based labels with SingleR

Markers and a dot plot are enough for broad cell classes. SingleR is an
optional second check.

```r
library(SingleR)
library(celldex)

ref <- celldex::MouseRNAseqData()
sort(unique(ref$label.main))           # which labels can it assign?

pred <- SingleR(test = GetAssayData(obj, layer = "data"),
                ref = ref, labels = ref$label.main)
obj$singler <- pred$pruned.labels      # NA = not a confident label

table(Idents(obj), obj$singler, useNA = "ifany")
plotScoreHeatmap(pred)
```

SingleR can only assign labels that are in the reference. A cell type missing
from the reference gets the closest available label.

## Save and record

```r
saveRDS(obj, "my_object.rds")
obj <- readRDS("my_object.rds")
sessionInfo()
```

Record every choice and its source in your decision record as you go.