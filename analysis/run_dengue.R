
library(HDF5Array)
library(SingleCellExperiment)
library(IsoCellChat)
library(zellkonverter)

file_path <- "data/GSE116672_only.h5"

# CREATE SCE OBJECT
expression_data <- HDF5Array(file_path, name = "data/expression")

transcripts_df <- read.csv("../../transcript_annotations.csv", stringsAsFactors = FALSE)
metadata_df    <- read.csv("../../refined_metadata.csv", stringsAsFactors = FALSE)


sce <- SingleCellExperiment(
  # 1. The Matrix: We name it 'counts' assuming it is raw data.
  # Transpose table
  assays = list(counts = t(expression_data)),

  # 2. The Columns:
  colData = metadata_df,

  # 3. The Rows:
  rowData = transcripts_df
)

rownames(sce) <- transcripts_df$transcript_id
colnames(sce) <- metadata_df$geo_accession


sample_list = c("1-026", "1-036", "1-013", "1-010", "1-008", "1-020", "3-013", "3-027", "3-018", "3-006")

for (patient_id in sample_list) {
  print(patient_id)
  sce_subset <- sce[, sce$patient.id == patient_id]
  sce_subset <- sce_subset[, sce_subset$cell_type_refined != "unknown"]

  results_list <- run_mega_cellchat_pipeline(
    sce = sce_subset,
    deeploc_file = "../ThesisProject/data/deeploc/deeploc_results.csv",
    species = "Human",
    assay_name = "counts",             # Your raw counts assay name
    gene_id_col = "gene_symbol",    # Your rowData symbol column
    cell_type_col = "cell_type_refined" # Your colData cell type column
  )

  results_list[["sce_baseline"]] <- sce_subset

  saveRDS(results_list,
          file=paste0("results/dengue_",patient_id,".rds"))
}




