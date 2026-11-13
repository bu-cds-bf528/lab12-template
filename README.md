## Single Cell Lab - QC Thresholds

In science, we often gloss over details that seem routine, like what we
filtered out before we began. Most papers report single cell QC in a single
sentence, and very few actually look at what was removed in the paper itself. 

In this lab, you'll take a published QC decision and work out what impact it
had on the dataset. The authors removed cells with fewer than 700 or more than
5,000 genes, or more than 20% mitochondrial reads. Your job is to decide whether
that was a defensible choice:

1. Make the strongest case for their filter.
2. Make the strongest case against it.
3. Give your rationale, including what it depends on and the evidence for each part.

To see what their filter removed, you'll start from a more lenient version of
the analysis.

## How to approach this lab

**Paper:** Joglekar A, et al. (2021). A spatially resolved brain region- and
cell type-specific isoform atlas of the postnatal mouse brain. *Nature
Communications* 12, 463. https://doi.org/10.1038/s41467-020-20343-5

Unlike the previous labs, this one has no step-by-step instructions. Start with
the paper: the methods give you the QC thresholds, the number of PCs, and the
clustering resolution, and the data availability statement tells you where the
data lives. Find the GEO accession there and download the feature matrices
(`barcodes`, `features`, and `matrix` files) for the two hippocampus replicates.
`Read10X()` expects the files in one folder per sample with exactly those
names, so you may need to rename them after downloading. Anything the paper
doesn't state, record as an assumption in your decision record.

## Deliverables

### 1. Code

- An R Markdown file that knits from top to bottom without errors
- Your analysis uses a lenient QC filter, so the cells the authors removed are
  still present, and you match the authors' number of PCs and clustering resolution

### 2. Decision record

- `decision_record.md`
- Mark which choices came from the paper's methods, and which you had to
  assume because the paper doesn't state them

### 3. Required figures and tables

At minimum:

- **QC distributions** for each sample, with the authors' cutoffs drawn on them
- **A UMAP of your clusters** with the cells the authors' filter would remove
  highlighted
- **A table** of how many cells the authors' filter removes from each cluster
- **Evidence for your annotation:** a marker plot, reference-based labels, or
  both, at the level of broad cell classes

### 4. Results

1. **What the filter removed.** Which cells, from which clusters or cell types,
   and how many.
2. **The case for the authors' thresholds.** At least three different kinds of
   justification (distribution, biology, precedent, consequences)

    - **Distribution:** Do the cutoffs fall at natural breaks in the data? Look
     for a dip, a shoulder, or a separate population in the QC distributions,
     and show where the authors' cutoffs sit relative to it.
   - **Biology:** Are the cutoffs reasonable for the cells this tissue should
     contain? Compare the gene counts and mitochondrial percentages of the cell
     types you identified with the cutoffs. Would healthy cells of each type
     pass?
   - **Precedent:** Are the cutoffs typical for this kind of data? Compare them
     with other common choices, such as standard tutorial defaults or a
     data-driven cutoff, and show how many cells each would remove from this
     dataset.
   - **Consequences:** What does the filter actually remove, and does losing
     those cells matter for the authors' goals? Show which clusters or cell
     types lose cells, and how many.

3. **The case against.** The strongest argument that the thresholds removed
   something they shouldn't have.
4. **Your judgement.** Was the filter defensible for the authors' goals? State
   what your answer depends on, and point to the evidence for each part.

### What you don't need

- Fine-grained annotation: broad cell classes are enough
- Integration of the two replicates
- A polished figure style: clear, labeled plots are enough