#' Clean Transcript IDs
#'
#' Removes version numbers and PAR regions from transcript IDs.
#' @param sce SingleCellExperiment object.
#' @return Cleaned SCE object.
#' @noRd
clean_transcript_ids_helper <- function(sce) {
  transcript_ids <- rownames(sce)
  cleaned_ids <- gsub("\\..*", "", transcript_ids)
  valid_transcripts <- !grepl("_PAR_Y$", transcript_ids)

  sce_cleaned <- sce[valid_transcripts, ]
  rownames(sce_cleaned) <- cleaned_ids[valid_transcripts]

  return(sce_cleaned)
}

#' Filter Protein Coding Transcripts
#'
#' Uses biomaRt to keep only protein-coding transcripts.
#' @param sce SingleCellExperiment object.
#' @param species "Human" or "Mouse".
#' @return Filtered SCE object.
#' @importFrom biomaRt useMart getBM
#' @noRd
filter_protein_coding_helper <- function(sce, species) {
  dataset_name <- if (species == "Human") "hsapiens_gene_ensembl" else "mmusculus_gene_ensembl"

  ensembl <- tryCatch({
    biomaRt::useMart("ensembl", dataset = dataset_name)
  }, error = function(e) {
    stop("Could not connect to BioMart. Check internet connection.")
  })

  transcript_ids <- rownames(sce)
  t_info <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_biotype"),
                           filters = "ensembl_transcript_id",
                           values = unique(transcript_ids),
                           mart = ensembl)

  pc_ids <- t_info$ensembl_transcript_id[t_info$transcript_biotype == "protein_coding"]
  sce[rownames(sce) %in% pc_ids, ]
}


#' Aggregate Transcripts to Gene Level
#'
#' Sums transcript counts to gene level using efficient rowsum.
#' @param sce SingleCellExperiment object.
#' @param assay_name Name of assay to aggregate.
#' @param gene_id_col Column name in rowData containing gene symbols.
#' @return Matrix of aggregated counts.
#' @importFrom SummarizedExperiment rowData assay
#' @noRd
aggregate_gene_helper <- function(sce, assay_name, gene_id_col) {
  gene_ids <- SummarizedExperiment::rowData(sce)[[gene_id_col]]

  # 1. Get the matrix
  mat <- SummarizedExperiment::assay(sce, assay_name)

  # 2. Handle Sparse Matrix (The source of your error)
  # rowsum works best on standard matrices.
  # If it's a sparse dgCMatrix, we convert to dense matrix first.
  # (Note: If your data is massive >100k cells, this step might need memory,
  #  but it is necessary for the standard 'rowsum' function).
  if (inherits(mat, "dgCMatrix") || inherits(mat, "sparseMatrix")) {
    mat <- as.matrix(mat)
  }

  # 3. Efficient Aggregation (Replaces slow 'aggregate' function)
  # rowsum sums the rows (transcripts) that belong to the same group (gene_ids)
  gene_expr_mat <- rowsum(mat, group = gene_ids)

  return(gene_expr_mat)
}
