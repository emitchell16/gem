# R/utilities/slim_sim_outputs.R
#
# One-off: rewrite finished simulation outputs with only the columns the post-sim code reads (+ zstd compression)
#   frees disk in <local temp>/gem_sims_<...>/sim=<label>/year=<y>/batch_<i>.parquet
#   results are unchanged: post-sim tables read only the kept columns; restart skip logic still sees every batch file
#
# use (in a fresh R session, with NO model run in progress):
#   source("<code dir>/R/setup_environment.R"); load_gem_packages(); source("<code dir>/R/configs.R")
#   paths <- get_paths(); source_gem_modules(paths, <project R folder>)   # defines hook_output_cols()
#   source("<code dir>/R/utilities/slim_sim_outputs.R")
#   slim_sim_outputs(dry_run = TRUE)    # sizes only, nothing changed
#   slim_sim_outputs(dry_run = FALSE)   # rewrite
#
#' kept columns: the same set the model saves (sim_output_cols(): core outcomes + hook_output_cols() of a project module)
#'   (any_of: absent columns are skipped). Not kept: RR scaling constants, random-number columns, working columns.
#' dc_calibration and psa folders are not touched (only sim=<label> folders directly under out_root)

slim_keep_cols <- function() {
  unique(c(
    # ids, grouping, weights
    "pid", "iso3", "w", "age", "female", "bmi", "bmi_class0",
    # outcomes (configs.R lists; cost_* columns are matched separately)
    all_outcomes_list, "obesity", "t2d", "cvd", "ckd", "eskd", "ac_death", "stroke_ever", "stroke_n", "cvd_event_n",
    "stroke", "cvd_event",
    # intervention columns (empty without a project module)
    hook_output_cols()))
}

slim_sim_outputs <- function(out_root = NULL, dry_run = TRUE, compression = "zstd") {
  # newest full-run folder if not given
  if (is.null(out_root)) {
    cand <- list.dirs(paths$local_temp, recursive = FALSE)
    cand <- cand[grepl("^gem_sims_", basename(cand))]
    if (!length(cand)) stop("no gem_sims_* folder in ", paths$local_temp)
    out_root <- cand[which.max(file.info(cand)$mtime)]
  }
  cat("Folder:", out_root, "\n")
  # guard: a model run writes temp_simulation_* folders; recent changes = run in progress
  tmp <- list.dirs(paths$local_temp, recursive = FALSE)
  tmp <- tmp[grepl("^temp_simulation_", basename(tmp))]
  if (length(tmp) && any(difftime(Sys.time(), file.info(tmp)$mtime, units = "mins") < 15))
    stop("a temp_simulation_* folder changed in the last 15 min: stop the model run first")
  if (!arrow::codec_is_available(compression)) { message(compression, " not available; using snappy"); compression <- "snappy" }

  keep <- slim_keep_cols()
  sims <- list.dirs(out_root, recursive = FALSE)
  sims <- sims[grepl("^sim=", basename(sims))]
  files <- list.files(sims, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)
  gb <- function(x) round(sum(file.size(x), na.rm = TRUE) / 1e9, 2)
  cat(length(sims), "sim folders,", length(files), "files,", gb(files), "GB\n")
  if (dry_run) {
    print(sort(sapply(sims, function(s) gb(list.files(s, recursive = TRUE, full.names = TRUE))), decreasing = TRUE))
    cat("dry run: nothing changed. Run with dry_run = FALSE to rewrite.\n")
    return(invisible(NULL))
  }

  before <- gb(files); n_done <- 0L; n_skip <- 0L
  for (k in seq_along(files)) {
    f <- files[k]
    d <- as.data.table(read_parquet(f, mmap = FALSE))                 # mmap = FALSE: no Windows file lock
    cols <- intersect(names(d), c(keep, grep("^cost_", names(d), value = TRUE)))
    if (length(cols) == ncol(d)) { n_skip <- n_skip + 1L; next }      # already slim (rerun safe)
    n0 <- nrow(d); w0 <- sum(d$w)
    d <- d[, ..cols]
    tmpf <- paste0(f, ".tmp")
    write_parquet(d, tmpf, compression = compression)
    chk <- as.data.table(read_parquet(tmpf, col_select = "w", mmap = FALSE))
    if (nrow(chk) != n0 || abs(sum(chk$w) - w0) > 1e-6 * max(1, w0)) {
      unlink(tmpf); stop("check failed (rows or weights differ): ", f, " - original kept")
    }
    rm(d, chk); gc(verbose = FALSE)
    if (unlink(f) != 0 || file.exists(f)) { unlink(tmpf); stop("could not delete original (file in use?): ", f) }
    if (!file.rename(tmpf, f)) stop("rename failed; slim copy left at: ", tmpf)
    n_done <- n_done + 1L
    if (k %% 50 == 0) cat(k, "of", length(files), "files\n")
  }
  after <- gb(list.files(sims, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE))
  cat("Rewritten:", n_done, "| already slim:", n_skip, "| size:", before, "->", after, "GB\n")
  invisible(c(before_gb = before, after_gb = after))
}
