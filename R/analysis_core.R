#' Run Standard CellChat (Gene Level)
#' @importFrom CellChat createCellChat subsetData identifyOverExpressedGenes identifyOverExpressedInteractions computeCommunProb computeCommunProbPathway aggregateNet normalizeData
#' @noRd
run_standard_cc <- function(expr_data, meta, cell_col, species) {

  # --- 1. SAFETY: Clean Metadata & Factors ---
  if (!cell_col %in% colnames(meta)) stop(paste("Column", cell_col, "not found in metadata"))

  meta[[cell_col]] <- as.factor(as.character(meta[[cell_col]]))
  meta[[cell_col]] <- droplevels(meta[[cell_col]])

  common_cells <- intersect(colnames(expr_data), rownames(meta))

  if (length(common_cells) == 0) stop("No matching cells found between Expression Data and Metadata!")
  if (length(common_cells) < ncol(expr_data) || length(common_cells) < nrow(meta)) {
    message("... Aligning Matrix and Metadata ...")
    expr_data <- expr_data[, common_cells, drop = FALSE]
    meta <- meta[common_cells, , drop = FALSE]
  }

  # --- SILENCE WARNINGS ---
  # Force 'samples' to be a FACTOR
  if (!"samples" %in% colnames(meta)) {
    meta$samples <- factor("sample1")
  } else {
    meta$samples <- as.factor(meta$samples)
  }

  # --- 2. Normalization ---
  expr_mat <- as.matrix(expr_data)
  expr_mat[is.na(expr_mat)] <- 0
  norm_mat <- CellChat::normalizeData(as(expr_mat, "CsparseMatrix"))

  # --- 3. Create Object ---
  cc <- CellChat::createCellChat(norm_mat, meta = meta, group.by = cell_col)
  cc@DB <- if(species == "Human") CellChat::CellChatDB.human else CellChat::CellChatDB.mouse

  # --- 4. Analysis Loop ---
  cc <- CellChat::subsetData(cc)
  cc <- CellChat::identifyOverExpressedGenes(cc, do.fast = FALSE)

  cc@meta[[cell_col]] <- droplevels(cc@meta[[cell_col]])
  cc@idents <- factor(cc@meta[[cell_col]])

  cc <- CellChat::identifyOverExpressedInteractions(cc)

  tryCatch({
    cc <- CellChat::computeCommunProb(cc, type = "truncatedMean", trim = 0.05)
  }, error = function(e) {
    message("Warning: truncatedMean failed. Switching to standard 'mean'.")
    cc <<- CellChat::computeCommunProb(cc, type = "mean")
  })

  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)

  print("CellChat analysis complete!")

  interactions <- subsetCommunication(cc)
  interaction_counts <- as.data.frame(as.table(cc@net$count))
  colnames(interaction_counts) <- c("Sender", "Receiver", "Count")

  return(list(cellchat = cc, interactions = interactions, interaction_counts = interaction_counts))
}

#' Run Isoform CellChat (Domain Aware)
#' @importFrom SummarizedExperiment assay colData
#' @noRd
run_isoform_cc <- function(sce, db, cell_col) {

  data.input <- SummarizedExperiment::assay(sce, "X")

  # --- DIAGNOSTIC CHECK (Temporary) ---
  message("--- DEBUG: Checking Data vs Database Compatibility ---")
  data_ids <- rownames(data.input)
  db_ids <- rownames(db$geneInfo)

  message("Head of SCE Rownames: ", paste(head(data_ids, 3), collapse = ", "))
  message("Head of DB Keys:      ", paste(head(db_ids, 3), collapse = ", "))

  overlap_count <- length(intersect(data_ids, db_ids))
  message("Total Overlapping IDs: ", overlap_count)

  if (overlap_count == 0) {
    stop("CRITICAL ERROR: No overlap found between Data and Database! subsetData will remove all features.")
  }
  # ------------------------------------

  if (max(data.input, na.rm=TRUE) > 50) data.input <- log1p(data.input)

  meta <- as.data.frame(SummarizedExperiment::colData(sce))

  if (!"samples" %in% colnames(meta)) {
    meta$samples <- factor("sample1")
  } else {
    meta$samples <- as.factor(meta$samples)
  }

  cc <- CellChat::createCellChat(data.input, meta = meta, group.by = cell_col)
  cc@DB <- db

  cc <- CellChat::subsetData(cc)

  # Safety check post-subsetData
  if (nrow(cc@data.signaling) == 0) {
    stop("subsetData removed all rows. Check the DEBUG output above for ID mismatch.")
  }
  cc <- CellChat::identifyOverExpressedGenes(cc, do.fast = FALSE)

  cc@meta[[cell_col]] <- droplevels(cc@meta[[cell_col]])
  cc@idents <- factor(cc@meta[[cell_col]])

  cc <- CellChat::identifyOverExpressedInteractions(cc)
  cc <- CellChat::computeCommunProb(cc, type = "truncatedMean", trim = 0.05)
  cc <- CellChat::computeCommunProbPathway(cc)
  cc <- CellChat::aggregateNet(cc)

  print("CellChat analysis complete!")

  interactions <- subsetCommunication(cc)
  interaction_counts <- as.data.frame(as.table(cc@net$count))
  colnames(interaction_counts) <- c("Sender", "Receiver", "Count")

  return(list(cellchat = cc, interactions = interactions, interaction_counts = interaction_counts))
}
