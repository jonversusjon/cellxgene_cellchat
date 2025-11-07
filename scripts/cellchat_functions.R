# cellchat_functions.R
# Drop-in helpers for SCE/HDF5-backed CellChat runs (no Seurat).
# - SCE load via zellkonverter::readH5AD(use_hdf5=TRUE)
# - Donor-aware disease handling (auto / donor / cell)
# - Build CellChat from matrix + meta
# - YAML-driven DB subset, sources/targets, prob type, trim, thresh, etc.

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(DelayedArray)
  library(HDF5Array)
  library(CellChat)
  library(Matrix)
  library(dplyr)
})


# ------------------------------------------------------------------
# Core utilities
# ------------------------------------------------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b
.now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
now <- .now  # Alias for compatibility with run_cellchat_batch.R
log_msg <- function(msg) cat(sprintf("[%s] %s\n", .now(), msg))

# Single normalization function for all uses
.norm <- function(x) tolower(trimws(as.character(x)))

.coerce_db_subset <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x) && length(x) == 1L && grepl(",", x)) {
    return(trimws(strsplit(x, ",")[[1]]))
  }
  if (is.character(x)) return(x)
  as.character(x)
}


# ------------------------------------------------------------------
# 1) Load & prepare .h5ad -> list(data, meta, cell_type_col)
#    data: genes x cells (HDF5-backed DelayedMatrix), meta: data.frame
# ------------------------------------------------------------------
load_and_prepare_anndata_cfg <- function(h5ad_path, cfg) {
  if (!file.exists(h5ad_path)) stop(sprintf("h5ad not found: %s", h5ad_path))

  # ensure these are available
  suppressPackageStartupMessages({
    library(zellkonverter)
    library(SingleCellExperiment)
    library(DelayedArray)
    library(HDF5Array)
    library(DelayedMatrixStats)
  })

  sce <- zellkonverter::readH5AD(h5ad_path, use_hdf5 = TRUE)

  # Expression (Delayed/HDF5-backed)
  assay_name <- if ("X" %in% SummarizedExperiment::assayNames(sce))
    "X" else SummarizedExperiment::assayNames(sce)[1]
  X <- SummarizedExperiment::assay(sce, assay_name)

  # Gene symbols
  genes <- as.character(SummarizedExperiment::rowData(sce)[["feature_name"]])
  rownames(X) <- make.unique(genes)

  # Metadata
  meta <- as.data.frame(SummarizedExperiment::colData(sce))
  rownames(meta) <- colnames(sce)

  # Choose cell-type column
  cell_type_col <- cfg$metadata$cell_type_col
  if (is.null(cell_type_col) || !(cell_type_col %in% colnames(meta))) {
    alt <- cfg$metadata$alt_cell_type_cols %||% c("labels","celltype","cluster","cell_type")
    hit <- alt[alt %in% colnames(meta)]
    if (length(hit)) cell_type_col <- hit[1]
  }
  if (is.null(cell_type_col) || !(cell_type_col %in% colnames(meta))) {
    stop("No cell type column found: set metadata.cell_type_col or add to alt_cell_type_cols.")
  }

  # Lightweight QC WITHOUT forcing in-memory conversion
  if (isTRUE(cfg$qc$drop_empty_genes) || !is.null(cfg$qc$min_cells_per_gene)) {
    nz_by_gene <- DelayedMatrixStats::rowSums2(X != 0)
    if (isTRUE(cfg$qc$drop_empty_genes)) {
      keep <- nz_by_gene > 0
      if (any(!keep)) X <- X[keep, , drop = FALSE]
    }
    if (!is.null(cfg$qc$min_cells_per_gene)) {
      thr <- as.integer(cfg$qc$min_cells_per_gene)
      keep <- nz_by_gene >= thr
      if (any(!keep)) X <- X[keep, , drop = FALSE]
    }
  }

  # IMPORTANT: do NOT coerce to dgCMatrix here; leave HDF5-backed
  list(
    data = X,
    meta = meta,
    cell_type_col = cell_type_col,
    assay_name = assay_name
  )
}

# Safely convert to dgCMatrix *after* sampling/subsetting
finalize_matrix_for_cellchat <- function(m) {
  if (inherits(m, "dgCMatrix")) return(m)

  suppressPackageStartupMessages(library(DelayedMatrixStats))
  # Check NNZ safely without realizing dense
  nnz <- sum(m != 0)
  # guard against 32-bit index limit for dgCMatrix
  if (nnz > .Machine$integer.max - 1L) {
    stop(sprintf(
      "Subset still has too many nonzeros for dgCMatrix (~%s). Reduce sampling/global cap or tighten QC.",
      format(nnz, big.mark = ",", scientific = FALSE)
    ))
  }
  # Convert via sparseMatrix to keep it sparse
  m <- as(m, "sparseMatrix")
  if (!inherits(m, "dgCMatrix")) m <- as(m, "dgCMatrix")
  m
}

# ------------------------------------------------------------------
# 2) Donor-aware disease handling: adds meta$donor_context
#    disease_level: "auto" | "donor" | "cell"
# ------------------------------------------------------------------
derive_donor_context <- function(meta, cfg, disease_col, normal_value) {
  donor_col <- cfg$metadata$donor_col %||% "donor_id"
  disease_level <- cfg$metadata$disease_level %||% "auto"     # "auto" | "donor" | "cell"
  disease_threshold <- cfg$metadata$disease_threshold %||% 0.1 # fraction threshold

  if (!(donor_col %in% colnames(meta))) stop(sprintf("derive_donor_context: missing donor_col '%s'", donor_col))
  if (!(disease_col %in% colnames(meta))) stop(sprintf("derive_donor_context: missing disease_col '%s'", disease_col))

  meta[[donor_col]]   <- as.character(meta[[donor_col]])
  meta[[disease_col]] <- as.character(meta[[disease_col]])

  if (identical(disease_level, "donor")) {
    meta$donor_context <- meta[[disease_col]]
    disease_encoding <- "donor_level"
    mixed_donors <- 0L
    total_donors <- dplyr::n_distinct(meta[[donor_col]])
    return(list(meta=meta, disease_encoding=disease_encoding,
                mixed_donors=mixed_donors, total_donors=total_donors))
  }

  if (identical(disease_level, "cell")) {
    tmp <- meta |>
      dplyr::group_by(.data[[donor_col]]) |>
      dplyr::mutate(
        disease_fraction = mean(.data[[disease_col]] != normal_value),
        donor_context = ifelse(disease_fraction > disease_threshold, "disease", normal_value)
      ) |>
      dplyr::ungroup()
    tmp$disease_fraction <- NULL

    donor_disease_df <- meta |>
      dplyr::select(dplyr::all_of(c(donor_col, disease_col))) |>
      dplyr::distinct() |>
      dplyr::count(.data[[donor_col]], name = "n_unique_diseases")

    mixed_donors <- sum(donor_disease_df$n_unique_diseases > 1)
    total_donors <- nrow(donor_disease_df)
    return(list(meta=tmp, disease_encoding="cell_level_aggregated",
                mixed_donors=mixed_donors, total_donors=total_donors))
  }

  # auto
  donor_disease_df <- meta |>
    dplyr::select(dplyr::all_of(c(donor_col, disease_col))) |>
    dplyr::distinct() |>
    dplyr::count(.data[[donor_col]], .data[[disease_col]], name = "n") |>
    dplyr::count(.data[[donor_col]], name = "n_unique_diseases")

  mixed_donors <- sum(donor_disease_df$n_unique_diseases > 1)
  total_donors <- nrow(donor_disease_df)

  if (mixed_donors == 0) {
    meta$donor_context <- meta[[disease_col]]
    disease_encoding <- "donor_level"
  } else if (mixed_donors < 0.05 * total_donors) {
    flag_df <- meta |>
      dplyr::group_by(.data[[donor_col]]) |>
      dplyr::summarise(any_disease = any(.data[[disease_col]] != normal_value), .groups="drop")
    meta <- meta |>
      dplyr::left_join(flag_df, by = donor_col) |>
      dplyr::mutate(donor_context = ifelse(any_disease, "disease", normal_value)) |>
      dplyr::select(-any_disease)
    disease_encoding <- "mostly_donor"
  } else {
    meta <- meta |>
      dplyr::group_by(.data[[donor_col]]) |>
      dplyr::mutate(
        disease_fraction = mean(.data[[disease_col]] != normal_value),
        donor_context = ifelse(disease_fraction > disease_threshold, "disease", normal_value)
      ) |>
      dplyr::ungroup()
    meta$disease_fraction <- NULL
    disease_encoding <- "cell_level_threshold"
  }

  list(meta=meta, disease_encoding=disease_encoding,
       mixed_donors=mixed_donors, total_donors=total_donors)
}

# ------------------------------------------------------------------
# 3) Build CellChat core (no Seurat): matrix + meta
# ------------------------------------------------------------------
build_cellchat_core <- function(data.input,
                                meta,
                                sample_name,
                                species = "human",
                                cell_type_col,
                                database_subset = NULL,
                                do_ppi = FALSE,
                                ppi_species = species) {

  if (is.null(rownames(data.input)) || is.null(colnames(data.input)))
    stop("data.input must have rownames (genes) and colnames (cells).")
  if (is.null(rownames(meta)) || !all(colnames(data.input) %in% rownames(meta)))
    stop("meta rownames must include all cells in data.input.")
  if (!(cell_type_col %in% colnames(meta)))
    stop(sprintf("cell_type_col '%s' not in meta.", cell_type_col))

  # Extract cell type labels and ensure no unused factor levels
  cell_labels <- meta[colnames(data.input), cell_type_col, drop = TRUE]
  
  # If it's a factor, drop unused levels to avoid CellChat errors
  if (is.factor(cell_labels)) {
    cell_labels <- droplevels(cell_labels)
  }
  
  # Convert to character to avoid factor-related issues in CellChat
  cell_labels <- as.character(cell_labels)
  
  # Check unique groups - allow single group for self-communication analysis
  unique_groups <- unique(cell_labels)
  if (length(unique_groups) < 1) {
    stop("No cell groups found in the data")
  }
  
  # Note: CellChat can handle single-group analysis (e.g., neurons to neurons)
  # The communication analysis will examine interactions within that cell type

  # Create metadata for CellChat with labels and samples
  # Use donor_id as samples if available for replicate-aware analysis
  meta_cc <- data.frame(labels = cell_labels,
                        row.names = colnames(data.input), check.names = FALSE)
  
  # Add samples column (donor_id) if available for replicate-aware analysis
  if ("donor_id" %in% colnames(meta)) {
    donor_ids <- as.character(meta[colnames(data.input), "donor_id"])
    meta_cc$samples <- donor_ids
  } else {
    # Fallback: all cells from same sample
    meta_cc$samples <- "sample1"
  }

  # choose DB
  db_name <- if (tolower(species) == "mouse") "CellChatDB.mouse" else "CellChatDB.human"
  CellChatDB <- get(db_name, envir = asNamespace("CellChat"))
  if (!is.null(database_subset)) {
    CellChatDB <- subsetDB(CellChatDB, search = database_subset)
  }

  # create object from matrix (suppress verbose cell group listing)
  invisible(capture.output({
    cellchat <- createCellChat(object = data.input, meta = meta_cc, group.by = "labels")
  }))
  cellchat@DB <- CellChatDB
  
  # Drop unused factor levels from @idents to avoid CellChat errors
  if (is.factor(cellchat@idents)) {
    cellchat@idents <- droplevels(cellchat@idents)
  }

  # subset to expressed, over-expression tests
  invisible(capture.output({
    cellchat <- subsetData(cellchat)
  }))
  
  # Drop unused factor levels again after subsetData
  if (is.factor(cellchat@idents)) {
    cellchat@idents <- droplevels(cellchat@idents)
  }
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # optional PPI step
  if (isTRUE(do_ppi)) {
    # Placeholder: wire in your preferred PPI step here if you use it.
    # cellchat <- netEmbedding(cellchat, type = "functional")  # example
  }

  cellchat
}

# ------------------------------------------------------------------
# 4) Run comm-prob + filters + pathways + aggregation (+ centrality)
# ------------------------------------------------------------------
compute_comm <- function(cellchat,
                         type = c("triMean","truncatedMean"),
                         trim = 0.1,
                         min_cells_for_comm = 10,
                         use_ppi = FALSE,
                         compute_centrality = TRUE,
                         sources.use = NULL,
                         targets.use = NULL,
                         thresh = 0.05) {  # kept for backward compatibility but not used

  type <- match.arg(type)

  # Ensure no unused factor levels in @idents before computeCommunProb
  # This is critical to avoid "Please check unique(object@idents)" errors
  if (is.factor(cellchat@idents)) {
    cellchat@idents <- droplevels(cellchat@idents)
  }

  # computeCommunProb does NOT accept sources.use/targets.use in standard CellChat
  # These parameters should be used for subsetting before or filtering after
  # Note: 'thresh' parameter is no longer supported by current CellChat version
  # Kept as function parameter for backward compatibility but not passed to computeCommunProb
  cellchat <- computeCommunProb(
    cellchat,
    type = type,
    trim = trim
  )

  cellchat <- filterCommunication(cellchat, min.cells = min_cells_for_comm)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  if (isTRUE(compute_centrality)) {
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
  }
  cellchat
}

# ------------------------------------------------------------------
# 5) YAML-aware convenience wrapper used by your batch worker
# ------------------------------------------------------------------
run_cellchat_cfg <- function(data.input, meta, sample_name, cfg, cell_type_col) {
  species <- cfg$run$species %||% "human"

  db_subset <- NULL
  if (!is.null(cfg$database$subset)) db_subset <- .coerce_db_subset(cfg$database$subset)

  cellchat <- build_cellchat_core(
    data.input = data.input,
    meta = meta,
    sample_name = sample_name,
    species = species,
    cell_type_col = cell_type_col,
    database_subset = db_subset,
    do_ppi = isTRUE(cfg$comm_prob$use_ppi),
    ppi_species = cfg$comm_prob$ppi_species %||% species
  )

  type   <- cfg$comm_prob$type   %||% "triMean"
  trim   <- as.numeric(cfg$comm_prob$trim   %||% 0.1)
  thresh <- as.numeric(cfg$comm_prob$thresh %||% 0.05)
  min_cells_for_comm <- as.integer(cfg$comm_prob$min_cells_for_comm %||% 10)

  srcs <- if (isTRUE(cfg$comm_prob$sources$enable)) cfg$comm_prob$sources$cell_types else NULL
  tgts <- if (isTRUE(cfg$comm_prob$targets$enable)) cfg$comm_prob$targets$cell_types else NULL

  cellchat <- compute_comm(
    cellchat,
    type = type,
    trim = trim,
    min_cells_for_comm = min_cells_for_comm,
    use_ppi = isTRUE(cfg$comm_prob$use_ppi),
    compute_centrality = isTRUE(cfg$qc$compute_centrality %||% TRUE),
    sources.use = srcs,
    targets.use = tgts,
    thresh = thresh
  )

  cellchat
}

# ------------------------------------------------------------------
# 6) Save helpers
# ------------------------------------------------------------------
save_cellchat_results_qc <- function(cellchat, outdir, sample_name, cfg = NULL) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  rds_path <- file.path(outdir, paste0(sample_name, "_cellchat.rds"))
  saveRDS(cellchat, rds_path)

  # small net counts CSV (source-target edges > 0)
  if (!is.null(cellchat@net$count)) {
    net_counts <- cellchat@net$count
    df <- as.data.frame(as.table(net_counts))
    colnames(df) <- c("source","target","count")
    df <- df[df$count > 0, , drop = FALSE]
    utils::write.csv(df, file.path(outdir, paste0(sample_name, "_net_counts.csv")), row.names = FALSE)
  }
}

# ------------------------------------------------------------------
# Meta type assignment helpers
# ------------------------------------------------------------------
assign_meta_types <- function(cell_types, meta_types) {
  if (is.null(meta_types) || length(meta_types) == 0) return(as.character(cell_types))
  ct_norm <- .norm(cell_types)
  out <- rep(NA_character_, length(ct_norm))
  for (bucket in names(meta_types)) {
    set_norm <- .norm(meta_types[[bucket]])
    hit <- ct_norm %in% set_norm
    out[hit] <- bucket
  }
  # fallback to original label if no bucket matched
  out[is.na(out)] <- as.character(cell_types)[is.na(out)]
  out
}

# ------------------------------------------------------------------
# Cell type sampling to manage memory
# ------------------------------------------------------------------
apply_cell_type_sampling <- function(data, meta, cell_type_col, cfg, seed = 123) {
  # Extract sampling config
  sampling_cfg <- cfg$sampling
  if (is.null(sampling_cfg) || !isTRUE(sampling_cfg$enable)) {
    return(list(data = data, meta = meta, sampled = FALSE))
  }
  
  max_per_type <- as.integer(sampling_cfg$max_cells_per_type %||% Inf)
  global_cap <- as.integer(sampling_cfg$global_cap %||% 0)
  min_cells_to_sample <- as.integer(sampling_cfg$sample_only_if_ncells_gt %||% 0)
  
  if (max_per_type == Inf && global_cap == 0) {
    return(list(data = data, meta = meta, sampled = FALSE))
  }
  
  set.seed(seed)
  
  # Get cell types
  if (!(cell_type_col %in% colnames(meta))) {
    warning("Cell type column not found, skipping sampling")
    return(list(data = data, meta = meta, sampled = FALSE))
  }
  
  cell_types <- meta[[cell_type_col]]
  unique_types <- unique(cell_types)
  
  # Sample within each cell type
  keep_indices <- c()
  sampling_summary <- list()
  
  for (ct in unique_types) {
    ct_indices <- which(cell_types == ct)
    n_cells <- length(ct_indices)
    
    # Determine how many cells to keep
    if (n_cells <= min_cells_to_sample) {
      # Don't sample small cell types
      n_keep <- n_cells
      sampled_indices <- ct_indices
    } else if (n_cells <= max_per_type) {
      # Keep all if under the limit
      n_keep <- n_cells
      sampled_indices <- ct_indices
    } else {
      # Sample down to max_per_type
      n_keep <- max_per_type
      sampled_indices <- sample(ct_indices, n_keep, replace = FALSE)
    }
    
    keep_indices <- c(keep_indices, sampled_indices)
    sampling_summary[[as.character(ct)]] <- list(
      original = n_cells,
      kept = n_keep,
      sampled = n_keep < n_cells
    )
  }
  
  # Sort indices to maintain order
  keep_indices <- sort(keep_indices)
  
  # Apply global cap if specified
  if (global_cap > 0 && length(keep_indices) > global_cap) {
    keep_indices <- sample(keep_indices, global_cap, replace = FALSE)
    keep_indices <- sort(keep_indices)
  }
  
  # Subset data and metadata
  data_sampled <- data[, keep_indices, drop = FALSE]
  meta_sampled <- meta[keep_indices, , drop = FALSE]
  
  # Print summary
  n_original <- ncol(data)
  n_kept <- ncol(data_sampled)
  pct_kept <- round(100 * n_kept / n_original, 1)
  
  message(sprintf("Cell type sampling: %d/%d cells kept (%.1f%%)", 
                  n_kept, n_original, pct_kept))
  
  # Report per-type sampling
  sampled_types <- names(which(sapply(sampling_summary, function(x) x$sampled)))
  if (length(sampled_types) > 0) {
    message(sprintf("  Sampled cell types: %s", paste(sampled_types, collapse = ", ")))
    for (ct in head(sampled_types, 5)) {
      ss <- sampling_summary[[ct]]
      message(sprintf("    %s: %d -> %d cells", ct, ss$original, ss$kept))
    }
    if (length(sampled_types) > 5) {
      message(sprintf("    ... and %d more types", length(sampled_types) - 5))
    }
  }
  
  list(
    data = data_sampled, 
    meta = meta_sampled, 
    sampled = TRUE,
    summary = sampling_summary
  )
}

# Small helper to print overlaps per bucket
log_meta_overlap <- function(cell_types, meta_types, msgf=message) {
  if (is.null(meta_types) || length(meta_types) == 0) {
    msgf("Meta overlap: (no meta types loaded)")
    return(invisible(NULL))
  }
  ct_norm <- .norm(cell_types)
  counts <- vapply(names(meta_types), function(b) {
    sum(ct_norm %in% .norm(meta_types[[b]]))
  }, integer(1))
  top <- sort(counts, decreasing = TRUE)
  msgf(sprintf("Meta overlap counts (top 10): %s",
               paste(sprintf("%s=%d", names(top)[seq_len(min(10, length(top)))], 
                             top[seq_len(min(10, length(top)))]), collapse=", ")))
  if (sum(counts) == 0L) {
    uniq <- sort(unique(head(ct_norm, 40)))
    msgf(sprintf("No overlaps; first observed cell_type labels (normalized): %s",
                 paste(uniq, collapse=", ")))
  }
  invisible(counts)
}


# ------------------------------------------------------------------
# Fast meta/bucket lookup (using consolidated .norm function)
# ------------------------------------------------------------------

# Build once per subset
make_meta_lookup <- function(cell_types, meta_types) {
  ct_norm <- .norm(cell_types)

  # map raw label -> indices
  label_to_idx <- split(seq_along(ct_norm), ct_norm)

  # map bucket name -> indices (union of its member labels)
  bucket_to_idx <- list()
  if (!is.null(meta_types) && length(meta_types) > 0) {
    bucket_names_norm <- .norm(names(meta_types))
    for (i in seq_along(meta_types)) {
      bname <- bucket_names_norm[i]
      members <- .norm(meta_types[[i]])
      idx <- unlist(label_to_idx[names(label_to_idx) %in% members], use.names = FALSE)
      bucket_to_idx[[bname]] <- sort(unique(idx))
    }
  }

  list(
    ct_norm = ct_norm,
    label_to_idx = label_to_idx,         # named by normalized label
    bucket_to_idx = bucket_to_idx        # named by normalized bucket
  )
}

# Resolve a wanted vector (bucket names and/or raw labels) into indices
wanted_indices <- function(wanted, lookup) {
  if (is.null(wanted) || length(wanted) == 0) return(integer(0))

  w_norm <- .norm(wanted)
  out <- integer(0)
  # 1) bucket hits
  hit_b <- w_norm[w_norm %in% names(lookup$bucket_to_idx)]
  if (length(hit_b)) {
    out <- c(out, unlist(lookup$bucket_to_idx[hit_b], use.names = FALSE))
  }
  # 2) raw label hits
  hit_l <- w_norm[w_norm %in% names(lookup$label_to_idx)]
  if (length(hit_l)) {
    out <- c(out, unlist(lookup$label_to_idx[hit_l], use.names = FALSE))
  }
  sort(unique(out))
}