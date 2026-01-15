#' Run Complete Isoform CellChat Pipeline
#'
#' Executes the 4-stage analysis pipeline: Baseline, Protein-Coding, DeepLoc-Weighted, and Domain-Aware.
#'
#' @param sce Initial SingleCellExperiment object (Smart-Seq2 data).
#' @param deeploc_file Path to the DeepLoc CSV output file.
#' @param ppi_file Path to the Transcript-Specific PPI TSV. If NULL, looks for package default.
#' @param species Species name ("Human" or "Mouse"). Default "Human".
#' @param assay_name Name of the raw counts assay in SCE. Default "X".
#' @param gene_id_col Column name in rowData containing gene symbols. Default "gene_name".
#' @param cell_type_col Column name in colData containing cell types. Default "free_annotation".
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{sce_pc}: SCE object after protein-coding filtering.
#'   \item \code{sce_deeploc}: SCE object after DeepLoc weighting.
#'   \item \code{cc_baseline}: CellChat object (Baseline Gene Level).
#'   \item \code{cc_pc}: CellChat object (Protein Coding Gene Level).
#'   \item \code{cc_deeploc}: CellChat object (DeepLoc Adjusted).
#'   \item \code{cc_domain}: CellChat object (Domain Aware Isoform Level).
#' }
#' @export
run_mega_cellchat_pipeline <- function(sce,
                                       deeploc_file,
                                       ppi_file = NULL,
                                       species = "Human",
                                       assay_name = "X",
                                       gene_id_col = "gene_name",
                                       cell_type_col = "free_annotation") {

  message("=== Starting IsoCellChat Mega Pipeline ===")

  # Handle PPI File
  if (is.null(ppi_file)) {
    ppi_file <- system.file("extdata", "transcript_specific_interactions_2.tsv", package = "IsoCellChat")
    if (ppi_file == "") stop("PPI file not provided and default not found in package.")
  }

  # --- 1. Clean IDs ---
  message(">> Step 1: Cleaning Transcript IDs and Filtering to CellChatDB Genes...")
  sce_clean <- clean_transcript_ids_helper(sce)

  meta <- as.data.frame(SummarizedExperiment::colData(sce_clean))

  # --- 2. Baseline Analysis ---
  message(">> Step 2: Running Baseline CellChat (All Transcripts)...")
  gene_data_base <- aggregate_gene_helper(sce_clean, assay_name, gene_id_col)
  cc_baseline <- run_standard_cc(gene_data_base, meta, cell_type_col, species)

  # --- 3. Filter Protein Coding ---
  message(">> Step 3: Filtering for Protein Coding Transcripts...")
  sce_pc <- filter_protein_coding_helper(sce_clean, species)

  # --- 4. PC Analysis ---
  message(">> Step 4: Running PC Filtered CellChat...")
  gene_data_pc <- aggregate_gene_helper(sce_pc, assay_name, gene_id_col)
  cc_pc <- run_standard_cc(gene_data_pc, meta, cell_type_col, species)

  # --- 5. DeepLoc Weighting ---
  message(">> Step 5: Applying DeepLoc Weights...")
  sce_deeploc <- sce_pc
  sce_deeploc <- apply_deeploc_weights(sce_deeploc, deeploc_file, species, gene_id_col, assay_name)

  # --- 6. DeepLoc Analysis ---
  message(">> Step 6: Running DeepLoc Weighted CellChat...")
  gene_data_dl <- aggregate_gene_helper(sce_deeploc, assay_name, gene_id_col)
  cc_deeploc <- run_standard_cc(gene_data_dl, meta, cell_type_col, species)

  # --- 7. Domain Aware Analysis ---
  #message(">> Step 7: Running Domain Aware (Isoform) Analysis...")

  #isoform_db <- build_isoform_db(ppi_file, species)

  # Hybridize (Add Cofactors to the Weighted SCE)
  #sce_hybrid <- augment_cofactors(sce_deeploc, isoform_db, species, assay_name)

  #cc_domain <- run_isoform_cc(sce_hybrid, isoform_db, cell_type_col)

  message("=== Pipeline Complete ===")

  return(list(
    sce_pc = sce_pc,
    sce_deeploc = sce_deeploc,
    cc_baseline = cc_baseline,
    cc_pc = cc_pc,
    cc_deeploc = cc_deeploc #,
    #cc_domain = cc_domain
  ))
}
