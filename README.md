# IsoCellChat

**IsoCellChat** is an R package designed to perform isoform-aware cell-cell communication analysis. Unlike standard tools that aggregate expression to the gene level, IsoCellChat retains transcript-level resolution, applies protein localization filtering (DeepLoc), and uses a domain-aware database to infer specific isoform-isoform interactions.

## 📦 Installation

You can install the development version of IsoCellChat from GitHub:

```r
# install.packages("devtools")
devtools::install_github("DionysiosDimitreas/IsoCellChat")
```

## 📋 Data Requirements

To run the pipeline successfully, your input data must meet these specific requirements:

### 1. SingleCellExperiment (SCE) Object
Your input must be a `SingleCellExperiment` object containing:
* **Assay:** A raw counts matrix (typically named `"X"` or `"counts"`).
* **RowData:** A column containing **Gene Symbols** (e.g., `TGFB1`, `TNF`).
    * **Critical:** Do not point to Ensembl IDs (e.g., `ENSG000...`) here; CellChat requires standard symbols.
* **ColData:** A column containing **Cell Types** (e.g., `B cell`, `T cell`, `Monocyte`).

### 2. DeepLoc Results (CSV)
A CSV file containing subcellular localization predictions for your transcripts. It should include columns for:
* `Transcript_ID`
* `Localization`
* `Membrane` prediction scores
* `Extracellular` prediction scores

### 3. Transcript PPI Database (Optional)
A TSV file containing transcript-level protein-protein interactions. If not provided, the package uses its internal default database.

## 🚀 Usage

The package provides a "Mega Pipeline" function that runs the entire analysis workflow in a single command.

```r
library(IsoCellChat)
library(SingleCellExperiment)

# 1. Load your data
# sce <- readRDS("path/to/your_data.rds")

# 2. Run the Pipeline
results <- run_mega_cellchat_pipeline(
  sce = sce,
  deeploc_file = "path/to/deeploc_results.csv",
  species = "Human",             # "Human" or "Mouse"
  assay_name = "counts",         # Name of raw counts assay in your SCE
  gene_id_col = "gene_symbol",   # CRITICAL: rowData column with Gene Symbols
  cell_type_col = "annotation"   # colData column with Cell Types
)

# 3. Access Results
# The output is a list containing objects from every stage:

# Standard Gene-level CellChat (Baseline)
cc_baseline <- results$cc_baseline 

# Protein-Coding Filtered Gene-level
cc_pc <- results$cc_pc

# DeepLoc Weighted Gene-level
cc_deeploc <- results$cc_deeploc

# Final Isoform-Level Analysis (Domain Aware)
# Note: Full integration is in progress. See the "Domain-Aware Analysis" section below.
# cc_isoform <- results$cc_domain
```

## 🔬 Pipeline Workflow

The `run_mega_cellchat_pipeline` function performs the following steps automatically:

* **Filtering:** Subsets the dataset to keep only genes known to be involved in cell-cell communication (Ligands, Receptors, Cofactors) to optimize memory usage.
* **Cleaning:** Removes version numbers (e.g., `.1`) and PAR regions from transcript IDs.
* **Baseline Analysis:** Aggregates all transcripts to the gene level and runs standard CellChat.
* **Protein-Coding Filter:** Retains only protein-coding transcripts using Ensembl BioMart.
* **PC Analysis:** Runs CellChat on the protein-coding filtered dataset.
* **DeepLoc Weighting:** Down-weights transcripts predicted to be retained intracellularly (not secreted/membrane-bound).
* **DeepLoc Analysis:** Runs CellChat on the DeepLoc-weighted dataset.
* **Domain-Aware Analysis (In Progress):**
    * The full integration of domain-aware isoform inference is currently being developed for this package.
    * **Run Separately:** For now, you can run this specific analysis using our standalone repository here:
    * **[https://github.com/DionysiosDimitreas/Domain-aware-Cell-Cell-Communication][(https://github.com/DionysiosDimitreas/Domain-aware-Cell-Cell-Communication)]**
## 🛠️ Troubleshooting

### "CRITICAL ERROR: 0 Matches found"
This means the pipeline cannot match your SCE genes to the CellChat database.
* **Check `gene_id_col`:** Ensure you are pointing to the column with **Gene Symbols** (e.g., "TGFB1"), not Ensembl IDs or other identifiers.
* **Check Species:** Ensure `species = "Human"` or `"Mouse"` matches your data's capitalization (Human="TGFB1", Mouse="Tgfb1").

### "subsetData resulted in an empty matrix"
This usually happens if the `rownames` of your object were lost during processing. The latest version of the package includes robust fixes to preserve row names during the cofactor augmentation step.
