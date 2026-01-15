#' Build Isoform DB
#'
#' Parses the PPI TSV and constructs a CellChatDB object using rigorous complex expansion.
#' @importFrom dplyr filter distinct pull bind_rows select mutate all_of
#' @importFrom utils read.delim data
#' @importFrom stats na.omit
#' @noRd
build_isoform_db <- function(ppi_file, species) {

  # --- 1. Load Data ---
  if (!file.exists(ppi_file)) stop("PPI file not found.")
  my_isoform_ppi_df <- utils::read.delim(ppi_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

  # Load Base DB
  if (species == "Human") {
    utils::data("CellChatDB.human", package = "CellChat", envir = environment())
    base_db <- CellChatDB.human
  } else {
    utils::data("CellChatDB.mouse", package = "CellChat", envir = environment())
    base_db <- CellChatDB.mouse
  }

  # --- 2. Define Complex Logic ---
  complex_db <- base_db[["complex"]]
  all_complex_names <- rownames(complex_db)

  complex_ligands <- unique(base_db$interaction$ligand[base_db$interaction$ligand %in% all_complex_names])
  complex_receptors <- unique(base_db$interaction$receptor[base_db$interaction$receptor %in% all_complex_names])

  parse_complex <- function(name) {
    if (name %in% rownames(complex_db)) {
      subunits <- stats::na.omit(unlist(complex_db[name, ]))
      subunits[subunits != ""]
    } else character(0)
  }

  # --- 3. Filter Base Interactions ---
  genes_in_my_ppi <- unique(c(my_isoform_ppi_df$ligand_gene, my_isoform_ppi_df$receptor_gene))

  complex_subunits <- unique(c(complex_db$subunit_1, complex_db$subunit_2,
                               complex_db$subunit_3, complex_db$subunit_4, complex_db$subunit_5))
  complex_subunits_in_ppi <- complex_subunits[complex_subunits %in% genes_in_my_ppi]

  filtered_interactions <- dplyr::filter(base_db$interaction,
                                         (ligand %in% genes_in_my_ppi) |
                                           (receptor %in% genes_in_my_ppi) |
                                           (ligand %in% complex_ligands & ligand %in% complex_db$complex_name[complex_db$complex_name %in% complex_subunits_in_ppi]) |
                                           (receptor %in% complex_receptors & receptor %in% complex_db$complex_name[complex_db$complex_name %in% complex_subunits_in_ppi])
  )

  message(paste("Base interactions filtered to:", nrow(filtered_interactions)))

  # --- 4. Expansion Logic ---
  expand_cellchat_with_complex_handling <- function(cellchat_row) {

    ligand_symbol <- cellchat_row$ligand
    receptor_symbol <- cellchat_row$receptor
    ligand_is_complex <- ligand_symbol %in% complex_ligands
    receptor_is_complex <- receptor_symbol %in% complex_receptors

    get_isoforms_from_symbol <- function(gene_sym) {
      lig_sub <- dplyr::filter(my_isoform_ppi_df, ligand_gene == gene_sym)
      lig_ensts <- dplyr::pull(dplyr::distinct(lig_sub, ligand_transcript), ligand_transcript)
      rec_sub <- dplyr::filter(my_isoform_ppi_df, receptor_gene == gene_sym)
      rec_ensts <- dplyr::pull(dplyr::distinct(rec_sub, receptor_transcript), receptor_transcript)

      all_ensts <- unique(c(lig_ensts, rec_ensts))

      if (length(all_ensts) > 0) {
        return(data.frame(gene_symbol = gene_sym, ENST = all_ensts, ENSP = NA, stringsAsFactors = FALSE))
      } else {
        return(data.frame(gene_symbol = gene_sym, ENST = gene_sym, ENSP = NA, stringsAsFactors = FALSE))
      }
    }

    find_interacting_subunit_isoforms <- function(subunit_gene, partner_ensts) {
      sub_A <- dplyr::filter(my_isoform_ppi_df, ligand_gene == subunit_gene & receptor_transcript %in% partner_ensts)
      matches_A <- dplyr::pull(sub_A, ligand_transcript)
      sub_B <- dplyr::filter(my_isoform_ppi_df, receptor_gene == subunit_gene & ligand_transcript %in% partner_ensts)
      matches_B <- dplyr::pull(sub_B, receptor_transcript)
      return(unique(c(matches_A, matches_B)))
    }

    # CASE 1: Simple -> Simple
    if (!ligand_is_complex && !receptor_is_complex) {
      ligand_mapping <- get_isoforms_from_symbol(ligand_symbol)
      receptor_mapping <- get_isoforms_from_symbol(receptor_symbol)

      expanded_rows <- do.call(rbind, lapply(seq_len(nrow(ligand_mapping)), function(i) {
        ligand_enst <- ligand_mapping$ENST[i]
        ligand_ensp <- ligand_mapping$ENSP[i]

        do.call(rbind, lapply(seq_len(nrow(receptor_mapping)), function(j) {
          receptor_enst <- receptor_mapping$ENST[j]
          receptor_ensp <- receptor_mapping$ENSP[j]

          interaction_exists <- any(
            (my_isoform_ppi_df$ligand_transcript == ligand_enst & my_isoform_ppi_df$receptor_transcript == receptor_enst) |
              (my_isoform_ppi_df$receptor_transcript == ligand_enst & my_isoform_ppi_df$ligand_transcript == receptor_enst)
          )
          if (!interaction_exists) return(NULL)

          new_row <- cellchat_row
          new_row$ligand <- ligand_enst
          new_row$receptor <- receptor_enst
          new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
          new_row$interaction_name_2 <- paste(ligand_enst, "-", receptor_enst)
          new_row$ligand_gene_symbol <- ligand_symbol
          new_row$receptor_gene_symbol <- receptor_symbol
          new_row$ligand_ensp <- ligand_ensp
          new_row$receptor_ensp <- receptor_ensp
          new_row
        }))
      }))
      return(list(interactions = expanded_rows, complex_mapping = data.frame()))
    }

    # CASE 2: Simple -> Complex
    if (!ligand_is_complex && receptor_is_complex) {
      ligand_mapping <- get_isoforms_from_symbol(ligand_symbol)
      ligand_ensts <- ligand_mapping$ENST
      subunit_symbols <- parse_complex(receptor_symbol)

      subunit_interacting_isoforms <- list()
      for (subunit_symbol in subunit_symbols) {
        interacting <- find_interacting_subunit_isoforms(subunit_symbol, ligand_ensts)
        if (length(interacting) == 0) return(list(interactions = data.frame(), complex_mapping = data.frame()))
        subunit_interacting_isoforms[[subunit_symbol]] <- interacting
      }
      isoform_combinations <- expand.grid(subunit_interacting_isoforms, stringsAsFactors = FALSE)

      interactions <- do.call(rbind, lapply(1:nrow(isoform_combinations), function(i) {
        complex_isoforms <- as.character(isoform_combinations[i, ])
        complex_receptor_name <- paste(complex_isoforms, collapse = "_")
        do.call(rbind, lapply(1:nrow(ligand_mapping), function(j) {
          ligand_enst <- ligand_mapping$ENST[j]
          ligand_ensp <- ligand_mapping$ENSP[j]

          new_row <- cellchat_row
          new_row$ligand <- ligand_enst
          new_row$receptor <- complex_receptor_name
          new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
          new_row$interaction_name_2 <- paste(ligand_enst, "-", complex_receptor_name)
          new_row$ligand_gene_symbol <- ligand_symbol
          new_row$receptor_gene_symbol <- receptor_symbol
          new_row$ligand_ensp <- ligand_ensp
          new_row$receptor_ensp <- NA
          new_row
        }))
      }))

      complex_mapping <- do.call(rbind, lapply(1:nrow(isoform_combinations), function(i) {
        complex_isoforms <- as.character(isoform_combinations[i, ])
        row <- data.frame(complex_name = paste(complex_isoforms, collapse = "_"), stringsAsFactors = FALSE)
        for (k in 1:5) row[[paste0("subunit_", k)]] <- if(k <= length(complex_isoforms)) complex_isoforms[k] else NA
        row
      }))
      return(list(interactions = interactions, complex_mapping = complex_mapping))
    }

    # CASE 3: Complex -> Simple
    if (ligand_is_complex && !receptor_is_complex) {
      receptor_mapping <- get_isoforms_from_symbol(receptor_symbol)
      receptor_ensts <- receptor_mapping$ENST
      subunit_symbols <- parse_complex(ligand_symbol)

      subunit_interacting_isoforms <- list()
      for (subunit_symbol in subunit_symbols) {
        interacting <- find_interacting_subunit_isoforms(subunit_symbol, receptor_ensts)
        if (length(interacting) == 0) return(list(interactions = data.frame(), complex_mapping = data.frame()))
        subunit_interacting_isoforms[[subunit_symbol]] <- interacting
      }
      isoform_combinations <- expand.grid(subunit_interacting_isoforms, stringsAsFactors = FALSE)

      interactions <- do.call(rbind, lapply(1:nrow(isoform_combinations), function(i) {
        complex_isoforms <- as.character(isoform_combinations[i, ])
        complex_ligand_name <- paste(complex_isoforms, collapse = "_")
        do.call(rbind, lapply(1:nrow(receptor_mapping), function(j) {
          receptor_enst <- receptor_mapping$ENST[j]
          receptor_ensp <- receptor_mapping$ENSP[j]

          new_row <- cellchat_row
          new_row$ligand <- complex_ligand_name
          new_row$receptor <- receptor_enst
          new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
          new_row$interaction_name_2 <- paste(complex_ligand_name, "-", receptor_enst)
          new_row$ligand_gene_symbol <- ligand_symbol
          new_row$receptor_gene_symbol <- receptor_symbol
          new_row$ligand_ensp <- NA
          new_row$receptor_ensp <- receptor_ensp
          new_row
        }))
      }))

      complex_mapping <- do.call(rbind, lapply(1:nrow(isoform_combinations), function(i) {
        complex_isoforms <- as.character(isoform_combinations[i, ])
        row <- data.frame(complex_name = paste(complex_isoforms, collapse = "_"), stringsAsFactors = FALSE)
        for (k in 1:5) row[[paste0("subunit_", k)]] <- if(k <= length(complex_isoforms)) complex_isoforms[k] else NA
        row
      }))
      return(list(interactions = interactions, complex_mapping = complex_mapping))
    }

    # CASE 4: Complex -> Complex
    if (ligand_is_complex && receptor_is_complex) {
      ligand_subunit_symbols <- parse_complex(ligand_symbol)
      receptor_subunit_symbols <- parse_complex(receptor_symbol)

      all_lig_ensts <- unlist(lapply(ligand_subunit_symbols, function(s) get_isoforms_from_symbol(s)$ENST))
      all_rec_ensts <- unlist(lapply(receptor_subunit_symbols, function(s) get_isoforms_from_symbol(s)$ENST))

      lig_sub_iso <- list()
      for (s in ligand_subunit_symbols) {
        matches <- find_interacting_subunit_isoforms(s, all_rec_ensts)
        if (length(matches) == 0) return(list(interactions = data.frame(), complex_mapping = data.frame()))
        lig_sub_iso[[s]] <- matches
      }
      rec_sub_iso <- list()
      for (s in receptor_subunit_symbols) {
        matches <- find_interacting_subunit_isoforms(s, all_lig_ensts)
        if (length(matches) == 0) return(list(interactions = data.frame(), complex_mapping = data.frame()))
        rec_sub_iso[[s]] <- matches
      }

      lig_iso_combos <- expand.grid(lig_sub_iso, stringsAsFactors = FALSE)
      rec_iso_combos <- expand.grid(rec_sub_iso, stringsAsFactors = FALSE)

      interactions <- do.call(rbind, lapply(1:nrow(lig_iso_combos), function(i) {
        lig_cplx <- paste(as.character(lig_iso_combos[i, ]), collapse = "_")
        do.call(rbind, lapply(1:nrow(rec_iso_combos), function(j) {
          rec_cplx <- paste(as.character(rec_iso_combos[j, ]), collapse = "_")
          new_row <- cellchat_row
          new_row$ligand <- lig_cplx
          new_row$receptor <- rec_cplx
          new_row$interaction_name <- paste(ligand_symbol, receptor_symbol, sep = "_")
          new_row$interaction_name_2 <- paste(lig_cplx, "-", rec_cplx)
          new_row$ligand_gene_symbol <- ligand_symbol
          new_row$receptor_gene_symbol <- receptor_symbol
          new_row$ligand_ensp <- NA
          new_row$receptor_ensp <- NA
          new_row
        }))
      }))

      lig_map <- do.call(rbind, lapply(1:nrow(lig_iso_combos), function(i) {
        comp <- as.character(lig_iso_combos[i, ])
        row <- data.frame(complex_name = paste(comp, collapse = "_"), stringsAsFactors = FALSE)
        for (k in 1:5) row[[paste0("subunit_", k)]] <- if(k <= length(comp)) comp[k] else NA
        row
      }))
      rec_map <- do.call(rbind, lapply(1:nrow(rec_iso_combos), function(i) {
        comp <- as.character(rec_iso_combos[i, ])
        row <- data.frame(complex_name = paste(comp, collapse = "_"), stringsAsFactors = FALSE)
        for (k in 1:5) row[[paste0("subunit_", k)]] <- if(k <= length(comp)) comp[k] else NA
        row
      }))

      return(list(interactions = interactions, complex_mapping = rbind(lig_map, rec_map)))
    }
    return(NULL)
  }

  # --- 5. Execute Expansion ---
  message("Expanding interactions... This may take a moment.")
  interaction_list <- split(filtered_interactions, seq(nrow(filtered_interactions)))

  results_list <- lapply(interaction_list, function(x) {
    tryCatch(expand_cellchat_with_complex_handling(x), error = function(e) NULL)
  })

  isoform_interactions <- dplyr::bind_rows(lapply(results_list, `[[`, "interactions"))

  # --- CRITICAL: Enforce Column Order (Matches Your Original) ---
  # This ensures 'evidence' is at index 9, 'annotation' at 10, etc.
  target_order <- c(
    "interaction_name", "pathway_name", "ligand", "receptor", "agonist", "antagonist",
    "co_A_receptor", "co_I_receptor", "evidence", "annotation", "interaction_name_2",
    "is_neurotransmitter", "ligand.symbol", "ligand.family", "ligand.location",
    "ligand.keyword", "ligand.secreted_type", "ligand.transmembrane", "receptor.symbol",
    "receptor.family", "receptor.location", "receptor.keyword", "receptor.surfaceome_main",
    "receptor.surfaceome_sub", "receptor.adhesome", "receptor.secreted_type",
    "receptor.transmembrane", "version", "ligand_gene_symbol", "receptor_gene_symbol",
    "ligand_ensp", "receptor_ensp"
  )

  # Only select columns that exist to prevent errors, but enforce order for those that do
  isoform_interactions <- dplyr::select(isoform_interactions, dplyr::all_of(target_order))

  isoform_complexes <- dplyr::bind_rows(lapply(results_list, `[[`, "complex_mapping"))

  if(nrow(isoform_complexes) > 0) {
    isoform_complexes <- dplyr::distinct(isoform_complexes, complex_name, .keep_all = TRUE)
    rownames(isoform_complexes) <- isoform_complexes$complex_name
    isoform_complexes$complex_name <- NULL
  }

  # --- 6. Build GeneInfo (Transcripts) ---
  ligand_map <- data.frame(Symbol = my_isoform_ppi_df$ligand_transcript, gene.symbol = my_isoform_ppi_df$ligand_gene, stringsAsFactors = FALSE)
  receptor_map <- data.frame(Symbol = my_isoform_ppi_df$receptor_transcript, gene.symbol = my_isoform_ppi_df$receptor_gene, stringsAsFactors = FALSE)
  new_geneInfo <- unique(rbind(ligand_map, receptor_map))

  new_geneInfo <- new_geneInfo[!is.na(new_geneInfo$Symbol) & new_geneInfo$Symbol != "", ]
  rownames(new_geneInfo) <- new_geneInfo$Symbol

  # --- 7. Add Cofactors (Exact match to your script) ---
  cofactors <- unique(unlist(base_db$cofactor))
  cofactor_genes <- cofactors[!is.na(cofactors) & cofactors != ""]

  if(length(cofactor_genes) > 0) {
    cofactor_info <- data.frame(
      Symbol = cofactor_genes,
      gene.symbol = cofactor_genes,
      stringsAsFactors = FALSE
    )
    rownames(cofactor_info) <- cofactor_genes
    new_geneInfo <- rbind(new_geneInfo, cofactor_info)
  }

  # --- 8. Final Return ---
  myIsoformCellChatDB <- list(
    interaction = isoform_interactions,
    complex = isoform_complexes,
    cofactor = base_db$cofactor,
    geneInfo = new_geneInfo
  )
  return(myIsoformCellChatDB)
}

#' Augment SCE with Cofactors
#'
#' Adds gene-level aggregated rows for cofactors to the transcript-level SCE.
#' @importFrom SingleCellExperiment SingleCellExperiment sizeFactors
#' @importFrom SummarizedExperiment rowRanges<- rowData<- colData assays assay assayNames assayNames<- rowData
#' @importFrom biomaRt useMart getBM
#' @importFrom methods is as
#' @noRd
project_cofactors <- function(sce_pc, sce_pc_loc, db, species, assay_name = "X") {

  augment_cofactors <- function(sce_pc, sce_pc_loc, db, species = "Human", assay_name = "counts") {

    # Validate inputs
    if (!assay_name %in% assayNames(sce_pc)) {
      stop(paste("Assay", assay_name, "not found in sce_pc object"))
    }

    # Get the counts matrix (aggregation should be on raw counts)
    counts_matrix <- assay(sce_pc, assay_name)

    # Connect to the Ensembl BioMart
    message("Connecting to Ensembl BioMart...")
    ensembl_mart <- useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = "hsapiens_gene_ensembl",
      host = "https://useast.ensembl.org"
    )

    # Extract unique cofactor genes
    cofactor_genes <- db$cofactor %>%
      unlist() %>%
      unique()

    genes_to_pool <- cofactor_genes[!is.na(cofactor_genes) & cofactor_genes != ""]

    if (length(genes_to_pool) == 0) {
      warning("No valid cofactor genes found in db$cofactor")
      return(sce_pc_loc)
    }

    message(paste("Querying BioMart for", length(genes_to_pool), "cofactor genes..."))

    # Get transcript-to-gene mappings in one batch query
    tx_map <- getBM(
      attributes = c('ensembl_transcript_id', 'hgnc_symbol'),
      filters = 'hgnc_symbol',
      values = genes_to_pool,
      mart = ensembl_mart
    )

    message(paste("Found", nrow(tx_map), "transcript-to-gene mappings"))

    # Aggregate transcripts for each gene
    pooled_rows <- list()
    genes_found <- 0
    genes_skipped <- 0

    message("Aggregating transcript counts to gene level...")

    for (gene in genes_to_pool) {
      # Find transcripts for this gene from BioMart mapping
      tx_for_gene <- tx_map$ensembl_transcript_id[tx_map$hgnc_symbol == gene]

      # Find which transcripts exist in the SCE object
      tx_in_sce <- intersect(tx_for_gene, rownames(counts_matrix))

      if (length(tx_in_sce) > 0) {
        # Extract counts for these transcripts
        tx_counts <- counts_matrix[tx_in_sce, , drop = FALSE]

        # Sum across transcripts - handle both single and multiple transcripts
        if (length(tx_in_sce) == 1) {
          # Single transcript: extract as vector
          pooled_expression <- as.vector(tx_counts)
        } else {
          # Multiple transcripts: sum them
          pooled_expression <- colSums(tx_counts)
        }

        pooled_rows[[gene]] <- pooled_expression
        genes_found <- genes_found + 1

        message(paste("  ✓", gene, ":", length(tx_in_sce), "transcript(s) aggregated"))
      } else {
        genes_skipped <- genes_skipped + 1
        warning(paste("  ✗", gene, ": No transcripts found in SCE object"))
      }
    }

    message(paste("\nSummary: Found", genes_found, "genes, skipped", genes_skipped))

    # Return unchanged object if no genes were pooled
    if (length(pooled_rows) == 0) {
      warning("No cofactor genes were successfully pooled. Returning original object.")
      return(sce_pc_loc)
    }

    # Create matrix of pooled gene-level counts
    new_rows_matrix <- do.call(rbind, pooled_rows)

    # --- FIX STARTS HERE ---

    # 1. Check if the original matrix is a DelayedArray
    # We use methods::is() to check inheritance safely
    is_delayed <- methods::is(counts_matrix, "DelayedArray")

    # 2. Match matrix structure
    if (is_delayed) {
      # If original is DelayedArray, the new one MUST be too
      requireNamespace("DelayedArray")
      new_rows_matrix <- DelayedArray::DelayedArray(new_rows_matrix)
    } else if (methods::is(counts_matrix, "sparseMatrix")) {
      # If original is sparse (and not Delayed), make new one sparse
      new_rows_matrix <- as(new_rows_matrix, "sparseMatrix")
    }

    # Match matrix type (sparse vs dense)
    if (is(counts_matrix, "sparseMatrix")) {
      new_rows_matrix <- as(new_rows_matrix, "sparseMatrix")
    }

    # Create temporary SCE for new rows
    sce_new_rows <- SingleCellExperiment(
      assays = stats::setNames(list(new_rows_matrix), assay_name),
      colData = colData(sce_pc_loc)
    )

    # Combine original and new rows
    message(paste("Combining:", nrow(sce_pc_loc), "original +",
                  nrow(sce_new_rows), "new features"))
    print(sce_pc_loc)
    print(sce_new_rows)
    sce_combined <- rbind(sce_pc_loc, sce_new_rows)

    message(paste("✓ Final SCE object:", nrow(sce_combined), "features"))
    message("⚠ WARNING: Re-normalize the object before downstream analysis!")

    return(sce_combined)
  }

  sce_hybrid <- augment_cofactors(results_list$sce_pc, results_list$sce_deeploc, myIsoformCellChatDB)

  return(sce_hybrid)
}
