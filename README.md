# Protein2Tissue
<!-- badges -->
![R](https://img.shields.io/badge/R-%3E%3D4.3.0-276DC3?logo=r)
![License](https://img.shields.io/badge/license-Artistic--2.0-green)
![Status](https://img.shields.io/badge/status-in%20development-orange)


`Protein2Tissue` is an R package for tissue enrichment analysis of
protein and gene lists. Given a set of proteins or genes of interest and an
optional background set, the package identifies tissues that are significantly
over-represented using statistical tests.

Results can be explored through multiple visualization and analysis tools,
including:
- Tissue specificity overviews
- Heatmaps
- Fisher’s exact test–based enrichment
- Bar charts
- Gene–tissue network plots
- Gene Ontology (GO) enrichment analysis

The bundled defaulty gene expression reference data are derived from the 
Human Protein Atlas (HPA), which classifies genes into tissue-specific 
expression categories.

---

## Tissue-specific expression categories

Protein2Tissue uses the following HPA tissue specificity classes:

- Tissue enriched — expression significantly higher in one tissue compared to all others
- Group enriched — elevated expression in a small group of tissues
- Tissue enhanced — expression higher in one tissue than the average across all tissues
- Not detected — no expression detected in any tissue

---

## Installation

`Protein2Tissue` is currently under development. Install the latest version
from GitHub:

```r

if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}
remotes::install_github("ChrYang/Protein2Tissue")
```

From Bioconductor (coming soon):

```r
BiocManager::install("Protein2Tissue")
```

---

## Overview

The main workflow of this package consists of 7 steps:

```
Input gene/protein list
        ↓
tissueAnalysis()        — map to Human Protein Atlas database
        ↓
plotProteinType()       — overview of RNA tissue specificity classes
        ↓
heatmapEnrich()         — gene × tissue expression fingerprint
        ↓
fisherByTissue()        — Fisher's exact test per tissue
        ↓
fisherPlot()            — visualise enrichment results
        ↓
networkPlot()           — gene-tissue bipartite network
        ↓
goEnriched()            — GO enrichment for a tissue of interest
```

---

## Loading the package
```r
library(tidyverse)
library(Protein2Tissue)
```
---

## Input data

### Reference database

The package ships with an example tissue expression database, ee_tb, derived
from the Human Protein Atlas.

Key columns include:
- Gene
- Ensembl
- tissue
- signal
- Tissue.specificity
- Is_secreted

---

## Scenario 1: Assessing tissue-specific protein expression

This workflow demonstrates tissue enrichment analysis using example plasma
proteomics data.

### Example input and background
```r
data("plasma_data", package = "Protein2Tissue")

myList <- plasma_data %>%
  filter(Coefficient.Age > 0 & q.Age < 0.05) %>%
  pull(UniqueSymbol) %>%
  toupper()

myBackground <- toupper(plasma_data$UniqueSymbol)
```

---

## Step 1 — Tissue enrichment analysis

`tissueAnalysis()` maps your input and background lists against the reference 
database (here we are using the shipped dataset in the package derived from the
Human Protein Atlas database) and returns an `Enrichment` S4 object:

```r
result <- tissueAnalysis(
    input      = myList,
    background = myBackground,
    typeKey    = "Gene",
    typeKeyBg  = "Gene"
)
```

Enrichment results can be saved to disk as CSV or TSV files:

```r
writeOutput(result,type = "csv","~/Enrichment_Results.csv")
```


---

## Step 2 — Tissue specificity overview

Before running statistical tests, it is useful to inspect how your input
genes are distributed across tissue specificity classes:

```r
plotProteinType(result)
```


This bar chart shows how many of your input genes fall into each specificity 
class (tissue enriched, group enriched, tissue enhanced,
low specificity, etc.), giving a quick overview of the tissue-specific
character of your gene list.

---

## Step 3 — Expression heatmap

`heatmapEnrich()` shows the tissue expression signal for each gene in your
input as a heatmap:

```r
heatmapEnrich(result)
```
---

## Step 4 — Fisher’s exact test

`fisherByTissue()` performs Fisher's exact test for each tissue to identify
which tissues are significantly enriched in your input list relative to the
background:

```r
fisherResult <- fisherByTissue(result, padj = "BH")
```
---

## Step 5 — Visualization

`fisherPlot()` generates a horizontal grouped bar chart comparing input
and background proportions across tissues, annotated with p-values:

```r
fisherPlot(fisherResult)
```

Blue bars show the proportion of input proteins expressed in each tissue;
grey bars show the background proportion. P-values from Fisher's exact test
are annotated to the right of each tissue.

---

## Step 6 — Network visualization


`networkPlot()` creates a bipartite network connecting genes to their
associated tissues, coloured by signal category:

```r
networkPlot(result)
```

In the network:

- **Large circles** are tissue hubs, labelled with the tissue name
- **Small points** are genes
- **Edge colour** encodes the signal category (tissue enriched, group
  enriched, or tissue enhanced)

This is particularly useful for identifying groups of genes that share
tissue associations.

---

## Step 7 — GO enrichment analysis

`goEnriched()` performs Gene Ontology enrichment analysis on genes
showing tissue-specific expression in a tissue of interest:

```r
goResult <- goEnriched(result, type = "BP", tissueTest = "Liver")
```
---

## Custom tissue expression databases

Custom tissue expression databases can be constructed from gene-by-tissue
expression matrices using classify_genes().(see below)



## Scenario 2: Tissue-Specific Expression Analysis of Differentially Expressed Genes

Protein2Tissue can also be used to identify genes or proteins enriched in a
specified tissue from a list of differentially expressed genes (DEGs) or
differentially expressed proteins (DEPs) derived from RNA‑seq or tissue‑specific
proteomics experiments.

In this mode, enrichment is assessed relative to a specific tissue instead of
across all tissues.

Example:

```r
result2 <- tissueAnalysis2(
  input         = c("BRCA1", "EGFR", "ETDB", 
                    "CCDC7", "TSPY9", "CCDC148"), #List of DEG/DEPs
  typeKey       = "Gene",
  database      = NULL,      # NULL uses bundled HPA database
  inclusion     = "All",
  secretoryOnly = FALSE,
  tissueToTest  = "Testis"   # tissue name must match the database
)
```

The results can be written to disk by the `write.csv()`  function.

```r
write.csv(result2,"~/result2.csv",row.names=F)
```

This approach is useful for:
- Validating tissue identity of DEGs
- Verifying tissue specificity in spatial or bulk omics experiments
- Targeted tissue hypothesis testing

---

## Scenario 3: Customising Your Own Tissue Expression Database

Protein2Tissue supports constructing a custom tissue expression database from
user‑supplied gene‑by‑tissue average expression matrices. This allows analysis
using alternative expression atlases, custom cohorts, or non‑human organisms.

In this example, a GTEx gene-by-tissue TPM matrix bundled with the package is
used.

---

### Classifying Tissue Specificity from an Expression Matrix

Genes are classified into tissue specificity categories using 
`geneClassification()`.

The input expression matrix must have:
- Rows: genes or proteins
- Columns: tissues
- Values: average expression (e.g. TPM)

Example:
```r
data("tissue_tpm_gtex")

res <- geneClassification(
  expression          = tissue_tpm_gtex,
  typeKey             = "Ensembl",
  min_expression      = 1,
  secretory_analysis  = FALSE
)
```
The output contains gene-level tissue specificity annotations compatible with
Protein2Tissue downstream workflows.

---

### Using Predefined Tissue Groups

Related tissues can be grouped prior to classification, allowing tissue
specificity to be evaluated at the group level rather than the automatically
ranked tissue levels.

Example:
```r
data("tissue_map")

res_grouped <- classify_genes(
  expression          = tissue_tpm_gtex,
  tissue_groups       = tissue_map,
  typeKey             = "Gene",
  min_expression      = 1,
  secretory_analysis  = FALSE
)
```
---

### Downstream Use of Custom Databases

The resulting classification table can be directly supplied to tissueAnalysis()
as a custom tissue expression database:

```r
result_custom <- tissueAnalysis(
  input    = myList,
  database = res_grouped
)
```

This enables tissue enrichment testing, visualization, and GO enrichment
analysis using user-defined tissue reference data.


---

## Dependencies

| package | use |
|---|---|
| `dplyr` | data manipulation |
| `ggplot2` | visualisation |
| `ggraph` | network visualisation |
| `igraph` | graph construction |
| `tidyr` | data reshaping |
| `clusterProfiler` | GO enrichment |
| `AnnotationDbi` | human gene annotation |
| `org.Hs.eg.db` | human gene annotation |
| `methods` | S4 class system |

---

## Data Source and Citation

Tissue specificity data are derived from the Human Protein Atlas (HPA).

If you use Protein2Tissue, please cite:

Uhlén M et al. (2015). Tissue-based map of the human proteome.
Science, 347(6220). https://doi.org/10.1126/science.1260419

Uhlén M et al. (2017). A pathology atlas of the human cancer transcriptome.
Science, 357(6352). https://doi.org/10.1126/science.aan2507

---

## License

MIT License — see [LICENSE](LICENSE) for details.
