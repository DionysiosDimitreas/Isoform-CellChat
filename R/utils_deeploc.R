#' Generate FASTA for DeepLoc
#'
#' Creates a FASTA file of protein sequences from the SCE object.
#' @param sce SingleCellExperiment object.
#' @param species "Human" or "Mouse".
#' @param output_file Path to save the FASTA file.
#' @export
generate_deeploc_fasta <- function(sce, species = "Human", output_file = "proteins.fasta") {

  sce_clean <- clean_transcript_ids_helper(sce)
  sce_pc <- filter_protein_coding_helper(sce_clean, species)

  ids <- rownames(sce_pc)
  dataset <- if (species == "Human") "hsapiens_gene_ensembl" else "mmusculus_gene_ensembl"
  mart <- biomaRt::useMart("ensembl", dataset = dataset)

  chunks <- split(ids, ceiling(seq_along(ids) / 500))
  seqs <- list()

  message("Querying sequences...")
  for(i in seq_along(chunks)) {
    res <- try(biomaRt::getBM(attributes=c("ensembl_transcript_id", "peptide"),
                              filters="ensembl_transcript_id", values=chunks[[i]], mart=mart), silent=TRUE)
    if(!inherits(res, "try-error")) seqs[[i]] <- res
  }

  all_seqs <- do.call(rbind, seqs)
  all_seqs <- all_seqs[all_seqs$peptide != "" & !is.na(all_seqs$peptide), ]
  all_seqs$peptide <- gsub("\\*$", "", all_seqs$peptide)

  lines <- paste0(">", all_seqs$ensembl_transcript_id, "\n", all_seqs$peptide)
  writeLines(lines, output_file)
  message(paste("FASTA saved to", output_file))
}

#' Apply DeepLoc Weights to Matrix
#'
#' Internal function to adjust counts based on localization probability.
#' @importFrom dplyr mutate filter inner_join left_join rename case_when select
#' @importFrom stringr str_detect str_remove_all str_extract str_trim
#' @importFrom scater logNormCounts normalizeCounts
#' @importFrom utils data
#' @importFrom SummarizedExperiment assay assay<- rowData
#' @noRd
apply_deeploc_weights <- function(sce, deeploc_path, species, gene_id_col, assay_name) {

  # --- 1. Define Target Locations Logic (Exact Match to Your Script) ---
  get_target_locations <- function(sp) {
    if (sp == "Human") {
      utils::data("CellChatDB.human", package = "CellChat", envir = environment())
      DB <- CellChatDB.human
    } else {
      utils::data("CellChatDB.mouse", package = "CellChat", envir = environment())
      DB <- CellChatDB.mouse
    }

    df <- DB$geneInfo

    # Safety check if column exists
    if (!"Location" %in% colnames(df)) return(NULL)

    # Select columns (Base R to avoid dplyr select issues)
    df <- df[, c("Symbol", "Location")]

    # --- EXACT CLEANING LOGIC ---
    # 1. Remove prefix
    df$Location <- stringr::str_remove_all(df$Location, "SUBCELLULAR LOCATION: ")
    # 2. Remove isoform tags
    df$Location <- stringr::str_remove_all(df$Location, "\\[.*?\\]:\\s*")
    # 3. Extract first location
    df$Location <- stringr::str_extract(df$Location, "^[^\\{\\.;,]+")
    # 4. Trim whitespace
    df$Location <- stringr::str_trim(df$Location)

    # --- EXACT CATEGORIZATION LOGIC ---
    df <- dplyr::mutate(df, Location_Category = dplyr::case_when(

      # Map specific organelles
      stringr::str_detect(Location, "Golgi") ~ "Golgi.apparatus",
      stringr::str_detect(Location, "Endoplasmic") ~ "Endoplasmic.reticulum",
      stringr::str_detect(Location, "Sarcoplasmic") ~ "Endoplasmic.reticulum",
      stringr::str_detect(Location, "Microsome") ~ "Endoplasmic.reticulum",
      stringr::str_detect(Location, "endoplasmic reticulum") ~ "Endoplasmic.reticulum",
      stringr::str_detect(Location, "Mitochondrion") ~ "Mitochondrion",
      stringr::str_detect(Location, "Peroxisome") ~ "Peroxisome",
      stringr::str_detect(Location, "Lysosome") ~ "Lysosome.Vacuole",
      stringr::str_detect(Location, "Vacuole") ~ "Lysosome.Vacuole",
      stringr::str_detect(Location, "Endosome") ~ "Lysosome.Vacuole",
      stringr::str_detect(Location, "endosome") ~ "Lysosome.Vacuole",

      # Map Nucleus
      stringr::str_detect(Location, "Nucleus") ~ "Nucleus",

      # Map Cell Membrane
      stringr::str_detect(Location, "Cell membrane") ~ "Cell.membrane",
      stringr::str_detect(Location, "cell membrane") ~ "Cell.membrane",
      stringr::str_detect(Location, "Cell surface") ~ "Cell.membrane",
      stringr::str_detect(Location, "Cell junction") ~ "Cell.membrane",
      stringr::str_detect(Location, "Synap") ~ "Cell.membrane",
      stringr::str_detect(Location, "Postsynaptic density") ~ "Cell.membrane",
      stringr::str_detect(Location, "Presynapse") ~ "Cell.membrane",
      stringr::str_detect(Location, "Myelin membrane") ~ "Cell.membrane",
      stringr::str_detect(Location, "Membrane raft") ~ "Cell.membrane",
      stringr::str_detect(Location, "Cell Projection") ~ "Cell.membrane",
      stringr::str_detect(Location, "Cell projection") ~ "Cell.membrane",
      stringr::str_detect(Location, "Membrane") ~ "Cell.membrane",
      stringr::str_detect(Location, "membrane") ~ "Cell.membrane",

      # Map Cytoplasm
      stringr::str_detect(Location, "Cytoplasm") ~ "Cytoplasm",
      stringr::str_detect(Location, "Cytosol") ~ "Cytoplasm",
      stringr::str_detect(Location, "Perikaryon") ~ "Cytoplasm",
      stringr::str_detect(Location, "Midbody") ~ "Cytoplasm",
      stringr::str_detect(Location, "Lipid droplet") ~ "Cytoplasm",
      stringr::str_detect(Location, "Melanosome") ~ "Cytoplasm",
      stringr::str_detect(Location, "granule") ~ "Cytoplasm",
      stringr::str_detect(Location, "Vesicle") ~ "Cytoplasm",

      # Map Extracellular
      stringr::str_detect(Location, "Secreted") ~ "Extracellular",
      stringr::str_detect(Location, "Extracellular") ~ "Extracellular",
      stringr::str_detect(Location, "Plastid") ~ "Plastid",

      # Specific Gene Overrides
      stringr::str_detect(Symbol, "CHAT") ~ "Cell.membrane",
      stringr::str_detect(Symbol, "ASMT") ~ "Cytoplasm",
      stringr::str_detect(Symbol, "DDC") ~ "Cytoplasm",
      stringr::str_detect(Symbol, "GAD1") ~ "Cell.membrane",
      stringr::str_detect(Symbol, "GLYCAM1") ~ "Extracellular",
      stringr::str_detect(Symbol, "HDC") ~ "Cytoplasm",
      stringr::str_detect(Symbol, "HGF") ~ "Extracellular",
      stringr::str_detect(Symbol, "KLK") ~ "Extracellular",
      stringr::str_detect(Symbol, "PNMT") ~ "Cytoplasm",
      stringr::str_detect(Symbol, "THBS") ~ "Cytoplasm",
      stringr::str_detect(Symbol, "TPH") ~ "Cytoplasm",
      stringr::str_detect(Symbol, "IFNA1") ~ "Extracellular",
      stringr::str_detect(Symbol, "CCL3L") ~ "Extracellular",
      stringr::str_detect(Symbol, "GPR1") ~ "Cell.membrane",

      TRUE ~ NA_character_
    ))

    return(df)
  }

  # --- 2. Map and Merge ---
  targets <- get_target_locations(species)
  if (is.null(targets)) return(sce)

  t_map <- data.frame(ensembl_transcript_id = rownames(sce), Symbol = SummarizedExperiment::rowData(sce)[[gene_id_col]])
  dl_out <- read.csv(deeploc_path)

  merged <- dplyr::rename(dl_out, ensembl_transcript_id = Protein_ID)
  merged <- dplyr::inner_join(merged, t_map, by = "ensembl_transcript_id")
  merged <- dplyr::left_join(merged, targets, by = "Symbol")

  # --- 3. Extract Probability (Matrix Lookup) ---
  # This replaces the 'rowwise() %>% get()' logic for package safety,
  # but performs the exact same mathematical operation.

  target_cols <- merged$Location_Category
  valid_rows <- !is.na(target_cols) & (target_cols %in% colnames(merged))

  if (any(valid_rows)) {
    col_indices <- match(target_cols[valid_rows], colnames(merged))
    row_indices <- which(valid_rows)
    extracted_probs <- merged[cbind(row_indices, col_indices)]
    merged$Probability <- NA_real_
    merged$Probability[valid_rows] <- as.numeric(extracted_probs)
  } else {
    merged$Probability <- NA_real_
  }

  merged <- dplyr::filter(merged, !is.na(Probability))
  merged <- merged[, c("ensembl_transcript_id", "Probability")]
  colnames(merged) <- c("Transcript_ID", "Probability")

  # --- 4. Apply Weights ---
  probs <- setNames(merged$Probability, merged$Transcript_ID)
  common <- intersect(rownames(sce), names(probs))

  if (length(common) == 0) return(sce)

  raw_mat <- SummarizedExperiment::assay(sce, assay_name)
  if (inherits(raw_mat, "dgCMatrix")) raw_mat <- as.matrix(raw_mat)

  raw_mat[common, ] <- raw_mat[common, ] * probs[common]
  SummarizedExperiment::assay(sce, assay_name) <- raw_mat

  return(sce)
}
