#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(yaml)
  library(future)
  library(future.apply)
})

# ------------------------------------------------------------------
# Configuration loading
# ------------------------------------------------------------------

# --- Args ---
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript run_cellchat_batch.R <config_yaml> [n_cores_override]\n")
  quit(status = 1)
}
CONFIG_PATH <- normalizePath(args[1])
if (!file.exists(CONFIG_PATH)) stop(sprintf("Config file not found: %s", CONFIG_PATH))
cfg <- yaml::read_yaml(CONFIG_PATH)

# --- Locate functions file (next to config) ---
CONFIG_DIR <- dirname(CONFIG_PATH)
FUNCTIONS_FILE <- file.path(CONFIG_DIR, "cellchat_functions.R")
if (!file.exists(FUNCTIONS_FILE)) {
  stop(sprintf("cellchat_functions.R not found at: %s", FUNCTIONS_FILE))
}
source(FUNCTIONS_FILE)

# --- Core resolution: CLI > config > SLURM > default(1) ---
slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "0"))
cfg_cores   <- if (!is.null(cfg$run$n_cores)) as.integer(cfg$run$n_cores) else NA_integer_
cli_cores   <- if (length(args) >= 2) as.integer(args[2]) else NA_integer_
N_CORES <- if (!is.na(cli_cores)) {
  log_msg(sprintf("Using CLI override: %d cores", cli_cores)); cli_cores
} else if (!is.na(cfg_cores)) {
  log_msg(sprintf("Using config setting: %d cores", cfg_cores)); cfg_cores
} else if (slurm_cores > 0) {
  log_msg(sprintf("Using SLURM allocation: %d cores", slurm_cores)); slurm_cores
} else {
  log_msg("No core specification found, using 1 core"); 1
}

# ----- Output dir override from environment (if provided) -----
out_env <- Sys.getenv("CELLCHAT_OUTDIR")  # returns "" if unset
if (nzchar(out_env)) {
  log_msg(sprintf("Overriding output directory with: %s", out_env))
  cfg$run$outdir <- out_env
}
stopifnot(!is.null(cfg$run$outdir), nzchar(cfg$run$outdir))

# create the run dir early; other code can assume it exists
dir.create(cfg$run$outdir, recursive = TRUE, showWarnings = FALSE)
OUTDIR <- normalizePath(cfg$run$outdir, mustWork = FALSE)

# create logs subdirectory
LOGSDIR <- file.path(OUTDIR, "logs")
dir.create(LOGSDIR, recursive = TRUE, showWarnings = FALSE)

# ---- Meta types load (mandatory if analyses use meta buckets) ----
meta_file <- cfg$metadata$meta_types_file %||% NA_character_
if (is.na(meta_file)) {
  # fallback to file next to the YAML
  default_meta <- file.path(CONFIG_DIR, "meta_cell_types.R")
  if (file.exists(default_meta)) {
    meta_file <- default_meta
    log_msg(sprintf("Using default meta_types_file: %s", meta_file))
  } else {
    log_msg("No meta_types_file set; analyses using meta buckets will match 0 cells.")
  }
}

# ------------------------------------------------------------------
# File resolution
# ------------------------------------------------------------------
.is_abs_path <- function(p) grepl("^(/|~)|^[A-Za-z]:[/\\\\]", p)

.resolve_input_files <- function(cfg) {
  base_dir <- cfg$input$dir %||% getwd()
  files_from_cfg <- cfg$input$files
  pattern <- cfg$input$pattern
  limit <- as.integer(cfg$input$limit_files %||% 0)

  if (!is.null(files_from_cfg) && length(files_from_cfg) > 0) {
    paths <- vapply(files_from_cfg, function(f) {
      f <- as.character(f)
      if (.is_abs_path(f)) path.expand(f) else file.path(base_dir, f)
    }, FUN.VALUE = character(1))
    paths <- normalizePath(paths[file.exists(paths)], mustWork = FALSE)
  } else if (!is.null(pattern)) {
    if (!dir.exists(base_dir)) stop(sprintf("Input dir does not exist: %s", base_dir))
    paths <- list.files(base_dir, pattern = pattern, full.names = TRUE)
  } else {
    stop("Must specify either 'files' or 'pattern' in input section")
  }

  paths <- paths[file.exists(paths)]
  if (length(paths) == 0) stop("No valid .h5ad files found")
  if (limit > 0 && length(paths) > limit) paths <- head(paths, limit)
  paths
}

h5ad_files <- .resolve_input_files(cfg)
log_msg(sprintf("Found %d .h5ad file(s)", length(h5ad_files)))
for (p in h5ad_files) log_msg(sprintf("  - %s", p))
log_msg(sprintf("Output dir: %s", normalizePath(OUTDIR, mustWork = FALSE)))
log_msg(sprintf("Parallel cores: %d", N_CORES))

# Load meta types from the designated file
meta_types <- NULL
if (!is.na(meta_file) && file.exists(meta_file)) {
  log_msg(sprintf("Loading meta types from: %s", meta_file))
  env <- new.env(parent = baseenv())
  sys.source(meta_file, envir = env)
  if (exists("define_meta_cell_types", envir = env, inherits = FALSE)) {
    meta_types <- env$define_meta_cell_types()
    # Stabilize and sanitize
    meta_types <- lapply(meta_types, function(v) sort(unique(as.character(v))))
    log_msg(sprintf("Meta types loaded: %d categories (%s...)",
                    length(meta_types), paste(head(names(meta_types), 3), collapse = ", ")))
  } else {
    log_msg(sprintf("WARNING: %s does not define define_meta_cell_types()", meta_file))
  }
} else {
  log_msg("No meta_types_file found; analyses using meta buckets will match 0 cells.")
}

# Analyses block (keeps your YAML-driven setup)
analyses <- cfg$cell_communication$analyses %||% list()
if (length(analyses) > 0) {
  log_msg(sprintf("Configured %d analyses", length(analyses)))
  for (a in analyses) {
    log_msg(sprintf("  - %s: %s -> %s",
                    a$name,
                    paste(a$senders, collapse=","),
                    paste(a$receivers, collapse=",")))
  }
}

# Disease comparison knobs
disease_comparison_cfg <- cfg$cell_communication$disease_comparison
run_disease_comparison <- isTRUE(disease_comparison_cfg$enable)
disease_col <- cfg$metadata$disease_col %||% "disease"
normal_value <- disease_comparison_cfg$normal_value %||% "normal"
min_cells_per_condition <- as.integer(disease_comparison_cfg$min_cells_per_condition %||% 10)
min_cells_per_analysis  <- as.integer(disease_comparison_cfg$min_cells_per_analysis %||% 20)

# ------------------------------------------------------------------
# Parallel plan
# ------------------------------------------------------------------
options(future.globals.maxSize = 64 * 1024^3)
if (length(h5ad_files) == 1L) {
  plan("sequential")
  log_msg("Using sequential plan (1 file) so stdout appears live.")
} else {
  n_workers <- min(N_CORES, length(h5ad_files))
  log_msg(sprintf("Setting up multisession plan with %d workers for %d files", n_workers, length(h5ad_files)))
  plan("multisession", workers = n_workers)
  log_msg(sprintf("Future plan configured: %s", class(plan())[1]))
  log_msg("About to call future_lapply - workers will spawn on first task...")
  
  # Test that future is working
  log_msg("Testing future workers with simple task...")
  test_result <- future_lapply(1:n_workers, function(i) {
    paste0("Worker ", i, " initialized on PID ", Sys.getpid())
  }, future.seed = TRUE)
  for (msg in test_result) log_msg(paste0("  ", msg))
  log_msg("Future workers test completed successfully!")
}

# ------------------------------------------------------------------
# Worker: SCE/HDF5-backed flow (NO Seurat)
# ------------------------------------------------------------------
log_msg("\n=== STARTING PARALLEL PROCESSING ===")
log_msg(sprintf("Processing %d files with %d worker(s)...", length(h5ad_files), if(length(h5ad_files) == 1L) 1 else n_workers))
log_msg("Note: With multisession, worker stdout is buffered until completion.")
log_msg("      Check system processes (ps/top) to see workers actively running.")
log_msg(sprintf("      Expected completion time: ~30-60 min per file (varies by size)"))
flush.console()

# Create a shared progress tracking file in logs/
progress_file <- file.path(LOGSDIR, "progress.txt")
cat(sprintf("[%s] BATCH_START: Processing %d files with %d worker(s)\n", 
            Sys.time(), length(h5ad_files), if(length(h5ad_files) == 1L) 1 else n_workers), 
    file = progress_file, append = FALSE)

results <- future_lapply(
  h5ad_files,
  future.seed = TRUE,
  future.globals = list(
    cfg = cfg,
    CONFIG_PATH = CONFIG_PATH,
    OUTDIR = OUTDIR,
    LOGSDIR = LOGSDIR,
    progress_file = progress_file,
    analyses = analyses,
    meta_types = meta_types,
    disease_col = disease_col,
    normal_value = normal_value,
    run_disease_comparison = run_disease_comparison,
    min_cells_per_condition = min_cells_per_condition,
    min_cells_per_analysis = min_cells_per_analysis,
    FUNCTIONS_FILE = FUNCTIONS_FILE
  ),
  FUN = function(h5) {
    # Helper function to log progress to shared file
    log_progress <- function(event, details = "") {
      timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      file_base <- sub("\\.h5ad$", "", basename(h5))
      msg_line <- sprintf("[%s] %s | %s | %s\n", timestamp, file_base, event, details)
      cat(msg_line, file = progress_file, append = TRUE)
    }
    
    # Log file start
    log_progress("FILE_START", sprintf("PID=%d", Sys.getpid()))
    
    suppressPackageStartupMessages({
      library(yaml)
      library(zellkonverter)
      library(SingleCellExperiment)
      library(DelayedArray)
      library(HDF5Array)
      library(CellChat)
      library(Matrix)
      library(dplyr)
    })
    source(FUNCTIONS_FILE)

    # Use log message helper with worker-specific prefix
    msg <- function(...) {
      cat(sprintf("[%s] [%s] %s\n", 
                  now(),  # Using the now() function from cellchat_functions.R
                  sub("\\.h5ad$", "", basename(h5)), 
                  paste0(..., collapse = "")))
    }

    # ---- Load meta types inside the worker (do not rely on captured globals) ----
    CONFIG_DIR <- dirname(CONFIG_PATH)  # derive from passed config path
    meta_file <- cfg$metadata$meta_types_file %||% file.path(CONFIG_DIR, "meta_cell_types.R")
    meta_types <- list()

    if (file.exists(meta_file)) {
      env <- new.env(parent = baseenv())
      sys.source(meta_file, envir = env)
      if (!exists("define_meta_cell_types", envir = env, inherits = FALSE)) {
        stop(sprintf("Meta file %s must define define_meta_cell_types()", meta_file))
      }
      meta_types <- env$define_meta_cell_types()
      # stabilize & sanitize
      meta_types <- lapply(meta_types, function(v) sort(unique(as.character(v))))
      msg(sprintf("WORKER: loaded meta types: %d buckets (%s…)",
                  length(meta_types), paste0(head(names(meta_types), 3), collapse = ", ")))
    } else {
      msg("WORKER: No meta_types_file; analyses using meta buckets will match 0 cells.")
    }

    # --- Core: Load and process data from h5ad ---
    msg("Loading .h5ad via zellkonverter (HDF5-backed)…")
    log_progress("DATA_LOAD_START", "")
    prepared <- load_and_prepare_anndata_cfg(h5, cfg)
    msg(sprintf("Loaded %d genes x %d cells", nrow(prepared$data), ncol(prepared$data)))
    log_progress("DATA_LOAD_COMPLETE", sprintf("%d genes x %d cells", nrow(prepared$data), ncol(prepared$data)))

    # --- Optional cell type sampling to manage memory ---
    # Apply sampling to manage memory for large datasets
    sampled_result <- apply_cell_type_sampling(
      data = prepared$data,
      meta = prepared$meta,
      cell_type_col = prepared$cell_type_col,
      cfg = cfg,
      seed = 123
    )
    prepared$data <- sampled_result$data
    prepared$meta <- sampled_result$meta
    was_sampled <- sampled_result$sampled

    # --- Donor-aware disease context (depends on loaded meta) ---
    disease_encoding <- "none"; mixed_donors <- 0L; total_donors <- 0L
    if (disease_col %in% colnames(prepared$meta)) {
      msg("Setting up donor-aware disease context…")
      context_res <- derive_donor_context(prepared$meta, cfg, disease_col, normal_value)
      prepared$meta <- context_res$meta
      disease_encoding <- context_res$disease_encoding
      mixed_donors <- context_res$mixed_donors
      total_donors <- context_res$total_donors
      msg(sprintf("  Disease encoding: %s, mixed_donors=%d, total_donors=%d",
                  disease_encoding, mixed_donors, total_donors))
    }

    # --- Define subset(s) for analyses ---
    subsets <- list(all = list(
      data   = prepared$data,
      meta   = prepared$meta,
      suffix = "_all",
      sender_filter = NULL  # no filtering for "all" subset
    ))

    if (run_disease_comparison && "donor_context" %in% colnames(prepared$meta)) {
      # add normal / disease subsets
      normal_mask  <- (prepared$meta$donor_context == normal_value)
      disease_mask <- (prepared$meta$donor_context != normal_value)
      
      msg(sprintf("Disease split: %d normal, %d disease cells", sum(normal_mask), sum(disease_mask)))
      
      # For normal subset: use only normal cells for both senders and receivers
      if (sum(normal_mask) > 0) {
        subsets$normal <- list(
          data   = prepared$data[, normal_mask, drop = FALSE],
          meta   = prepared$meta[normal_mask, , drop = FALSE],
          suffix = "_normal",
          sender_filter = NULL  # no additional filtering needed
        )
      }
      # For disease subset: use ALL cells but filter senders to disease only
      if (sum(disease_mask) > 0) {
        subsets$disease <- list(
          data   = prepared$data,  # ALL cells (disease + normal)
          meta   = prepared$meta,   # ALL cell metadata
          suffix = "_disease",
          sender_filter = which(disease_mask)  # indices of disease cells for sender filtering
        )
      }
    }

    # prepared$data : genes x cells (DelayedMatrix / HDF5Array)
    # prepared$meta : data.frame of cell metadata
    # prepared$cell_type_col : name of the cell-type column
    msg(sprintf("Prepared: %d genes x %d cells; cell_type_col='%s'",
                nrow(prepared$data), ncol(prepared$data), prepared$cell_type_col))

    # --- Run configured analyses per subset ---
    out_summaries <- list()
    for (sn in names(subsets)) {
      subd <- subsets[[sn]]
      cell_types <- subd$meta[[prepared$cell_type_col]]
      
      # === Build lookup once per subset (FAST) and check overlaps ===
      lookup <- make_meta_lookup(cell_types, meta_types)
      
      # Check overlap and assign meta types once (fixed duplicate code)
      if (!is.null(meta_types) && length(meta_types) > 0) {
        msg(sprintf("Checking meta overlap for subset='%s'…", sn))
        log_meta_overlap(cell_types, meta_types, msgf = function(x) msg(x))
        # Attach meta_type column once
        subd$meta$meta_type <- assign_meta_types(as.character(cell_types), meta_types)
      }

      for (a in analyses) {
        a_name <- a$name
        msg(sprintf("\n--- %s subset | Analysis: %s ---", sn, a_name))
        log_progress("ANALYSIS_START", sprintf("subset=%s analysis=%s", sn, a_name))

        # turn senders/receivers into index vectors quickly
        send_idx <- wanted_indices(a$senders, lookup)
        recv_idx <- wanted_indices(a$receivers, lookup)
        
        # For disease subset: filter senders to only disease cells
        if (!is.null(subd$sender_filter)) {
          n_before_filter <- length(send_idx)
          send_idx <- intersect(send_idx, subd$sender_filter)
          n_after_filter <- length(send_idx)
          if (n_before_filter > 0 && n_after_filter == 0) {
            msg(sprintf("  Applied disease sender filter: %d senders found in meta-type '%s', but 0 are from disease donors", 
                        n_before_filter, paste(a$senders, collapse=",")))
          } else {
            msg(sprintf("  Applied disease sender filter: %d/%d senders remain from disease donors", 
                        n_after_filter, n_before_filter))
          }
        }
        
        keep_idx <- sort(unique(c(send_idx, recv_idx)))

        n_send <- length(send_idx); n_recv <- length(recv_idx); n_keep <- length(keep_idx)
        msg(sprintf("  Senders=%d, Receivers=%d, Total used=%d", n_send, n_recv, n_keep))

        if (n_send < min_cells_per_condition || n_recv < min_cells_per_condition || n_keep < min_cells_per_analysis) {
          msg("  SKIPPED: insufficient cells for this analysis")
          log_progress("ANALYSIS_SKIPPED", sprintf("subset=%s analysis=%s reason=insufficient_cells senders=%d receivers=%d", 
                                                    sn, a_name, n_send, n_recv))
          out_summaries[[paste0(sn, "_", a_name)]] <- list(
            status = "skipped", reason = "insufficient_cells",
            n_senders = n_send, n_receivers = n_recv,
            subset = sn, analysis = a_name, file = basename(h5),
            disease_encoding = disease_encoding, mixed_donors = mixed_donors, total_donors = total_donors
          )
          next
        }

        # HDF5-backed slicing (genes x cells) by column indices
        data_use <- subd$data[, keep_idx, drop = FALSE]
        meta_use <- subd$meta[keep_idx, , drop = FALSE]
        
        # Ensure we have the same cell IDs after subsetting
        # Get the original cell IDs from the subset
        original_cell_ids <- rownames(subd$meta)[keep_idx]
        
        # Set consistent names on both data and meta
        if (!is.null(original_cell_ids)) {
          colnames(data_use) <- original_cell_ids
          rownames(meta_use) <- original_cell_ids
        }

        # ---- enforce name alignment strictly before CellChat core ----
        # colnames(data_use) must be identical to rownames(meta_use)
        if (is.null(colnames(data_use)) || is.null(rownames(meta_use))) {
          msg("WARNING: Cell IDs lost during subsetting, using indices as IDs")
          cell_ids <- paste0("cell_", seq_len(ncol(data_use)))
          colnames(data_use) <- cell_ids
          rownames(meta_use) <- cell_ids
        } else if (!identical(colnames(data_use), rownames(meta_use))) {
          msg(sprintf("WARNING: Name mismatch after subsetting. Data has %d cells, meta has %d rows", 
                      ncol(data_use), nrow(meta_use)))
          # Try to recover by using intersection
          inter <- intersect(colnames(data_use), rownames(meta_use))
          if (length(inter) > 0) {
            data_use <- data_use[, inter, drop = FALSE]
            meta_use <- meta_use[inter, , drop = FALSE]
            msg(sprintf("  Recovered %d matching cells", length(inter)))
          } else {
            stop("No matching cell IDs between data and metadata after subsetting")
          }
        }
        ##########################################################
        # DROP-IN REPLACEMENT START
        {
          # final assert (cheap)
          stopifnot(identical(colnames(data_use), rownames(meta_use)))

          # Convert to in-memory sparse *after* slicing (avoids SVT->dgC blowups)
          convert_after_subset <- function(m) {
            if (inherits(m, "dgCMatrix")) return(m)
            if (exists("finalize_matrix_for_cellchat", mode = "function")) {
              return(finalize_matrix_for_cellchat(m))
            }
            if (requireNamespace("DelayedMatrixStats", quietly = TRUE)) {
              nnz <- sum(m != 0)
              if (nnz > .Machine$integer.max - 1L) {
                stop(sprintf(
                  "Subset still has too many nonzeros (%s). Lower sampling/global cap or tighten QC.",
                  format(nnz, big.mark = ",", scientific = FALSE)
                ))
              }
            }
            m2 <- as(m, "sparseMatrix")
            if (!inherits(m2, "dgCMatrix")) m2 <- as(m2, "dgCMatrix")
            m2
          }
          data_use <- convert_after_subset(data_use)

          # tell CellChat exactly which groups are sources/targets (names, not indices)
          src_names <- unique(as.character(cell_types[send_idx]))
          tgt_names <- unique(as.character(cell_types[recv_idx]))
          all_groups <- union(src_names, tgt_names)

          sample_name <- paste0(sub("\\.h5ad$", "", basename(h5)), subd$suffix, "_", a_name)

          analysis_cfg <- cfg
          analysis_cfg$comm_prob$sources <- list(enable = TRUE, cell_types = src_names)
          analysis_cfg$comm_prob$targets <- list(enable = TRUE, cell_types = tgt_names)

          # Run with timing and warning capture
          wlist <- list()
          cellchat <- NULL
          withCallingHandlers({
            t0 <- Sys.time()
            cellchat <- tryCatch({
              run_cellchat_cfg(
                data.input   = data_use,
                meta         = meta_use,
                sample_name  = sample_name,
                cfg          = analysis_cfg,
                cell_type_col = prepared$cell_type_col
              )
            }, error = function(e) {
              msg(sprintf("  ERROR: %s", e$message))
              NULL
            })
            dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
            if (!is.null(cellchat)) {
              msg(sprintf("  ✓ %s finished in %.1fs (%d groups; %d senders, %d receivers)",
                          sample_name, dt, length(all_groups), n_send, n_recv))
              log_progress("ANALYSIS_CELLCHAT_SUCCESS", sprintf("subset=%s analysis=%s time=%.1fs groups=%d", 
                                                                  sn, a_name, dt, length(all_groups)))
            } else {
              log_progress("ANALYSIS_CELLCHAT_FAILED", sprintf("subset=%s analysis=%s", sn, a_name))
            }
          }, warning = function(w) {
            wlist[[length(wlist) + 1L]] <<- conditionMessage(w)
            invokeRestart("muffleWarning")
          })

          if (length(wlist)) {
            msg(sprintf("  ⚠ %d warnings (first 3): %s",
                        length(wlist), paste(utils::head(unique(wlist), 3L), collapse = " | ")))
          }
        }
        # DROP-IN REPLACEMENT END
        ##########################################################

        if (is.null(cellchat)) {
          out_summaries[[paste0(sn, "_", a_name)]] <- list(
            status = "failed", reason = "cellchat_error",
            n_senders = n_send, n_receivers = n_recv,
            subset = sn, analysis = a_name, file = basename(h5),
            disease_encoding = disease_encoding, mixed_donors = mixed_donors, total_donors = total_donors
          )
          next
        }

        save_cellchat_results_qc(cellchat, OUTDIR, sample_name, cfg)
        msg(sprintf("  DONE: saved results '%s'", sample_name))
        log_progress("ANALYSIS_COMPLETE", sprintf("subset=%s analysis=%s sample=%s interactions=%d", 
                                                    sn, a_name, sample_name, 
                                                    if (!is.null(cellchat@net$count)) sum(cellchat@net$count > 0) else 0))

        out_summaries[[paste0(sn, "_", a_name)]] <- list(
          status = "completed",
          n_senders = n_send, n_receivers = n_recv,
          n_interactions = if (!is.null(cellchat@net$count)) sum(cellchat@net$count > 0) else 0,
          subset = sn, analysis = a_name, file = basename(h5),
          disease_encoding = disease_encoding, mixed_donors = mixed_donors, total_donors = total_donors
        )

        rm(cellchat, data_use, meta_use); gc()
      }
    }

    rm(prepared); gc()
    log_progress("FILE_COMPLETE", sprintf("analyses_completed=%d", sum(sapply(out_summaries, function(x) x$status == "completed"))))
    out_summaries
  }
)

# ------------------------------------------------------------------
# Summary report
# ------------------------------------------------------------------
log_msg("\n=== GENERATING SUMMARY REPORT ===")
all_results <- do.call(c, results)
summary_df <- do.call(
  rbind,
  lapply(names(all_results), function(nm) {
    r <- all_results[[nm]]
    data.frame(
      analysis_name = nm,
      status        = r$status %||% NA,
      n_senders     = r$n_senders %||% NA,
      n_receivers   = r$n_receivers %||% NA,
      n_interactions= r$n_interactions %||% NA,
      subset        = r$subset %||% NA,
      analysis      = r$analysis %||% NA,
      file          = r$file %||% NA,
      stringsAsFactors = FALSE
    )
  })
)

write.csv(summary_df, file.path(OUTDIR, "analysis_summary.csv"), row.names = FALSE)

n_completed <- sum(summary_df$status == "completed", na.rm = TRUE)
n_skipped   <- sum(summary_df$status == "skipped",   na.rm = TRUE)
n_failed    <- sum(summary_df$status == "failed",    na.rm = TRUE)

log_msg(sprintf("\nFINAL SUMMARY:"))
log_msg(sprintf("  Total analyses attempted: %d", nrow(summary_df)))
log_msg(sprintf("  Completed: %d", n_completed))
log_msg(sprintf("  Skipped (insufficient cells): %d", n_skipped))
log_msg(sprintf("  Failed: %d", n_failed))

if (length(analyses) > 0) {
  log_msg("\nPer-analysis completion rates:")
  for (a in analyses) {
    n_runs <- sum(summary_df$analysis == a$name, na.rm = TRUE)
    n_succ <- sum(summary_df$analysis == a$name & summary_df$status == "completed", na.rm = TRUE)
    log_msg(sprintf("  %s: %d/%d successful", a$name, n_succ, n_runs))
  }
}

# Log batch completion
cat(sprintf("[%s] BATCH_COMPLETE: %d completed, %d skipped, %d failed\n", 
            Sys.time(), n_completed, n_skipped, n_failed), 
    file = progress_file, append = TRUE)

log_msg("All done.")