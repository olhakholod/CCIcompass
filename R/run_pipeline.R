#' Run the full CCIcompass pipeline from a Seurat object
#'
#' The main entry point for CCIcompass. Provide a Seurat object and the names
#' of your cell-type and donor columns; the function handles all LIANA calls,
#' donor-aware filtering, method intersection, and CIS scoring automatically,
#' returning a single ranked data frame ready for plotting.
#'
#' @details
#' The pipeline executes these steps in order:
#' \enumerate{
#'   \item Prepare Seurat identities (combine cell-type + donor labels if needed)
#'   \item Run LIANA with six CCI inference methods:
#'         connectome, logfc, natmi, sca, cellphonedb, cellchat
#'   \item Filter each method's output to same-donor source-target pairs
#'   \item Retain only interactions present in \code{min_donors} or more donors
#'   \item Intersect across all six methods (consensus interactions only)
#'   \item Compute CIS via Rank-Biased Precision (RBP) aggregation
#'   \item Final deduplication: one row per unique interaction
#' }
#'
#' \strong{Cell-type label format — two options:}
#' \itemize{
#'   \item \strong{Option A:} \code{idents_col} already embeds the donor ID
#'     after an underscore, e.g. \code{"Epi-Intestine_donor1"}.
#'     Set \code{donor_col = NULL} (default).
#'   \item \strong{Option B:} You have separate cell-type and donor columns.
#'     Set \code{donor_col} to the name of the donor column and the function
#'     creates the combined label automatically.
#' }
#'
#' @param seurat_obj A Seurat object. Required.
#' @param idents_col Character. Name of the \code{meta.data} column containing
#'   cell-type labels. Default \code{"celltype.combo"}.
#' @param donor_col Character or \code{NULL}. Name of the \code{meta.data}
#'   column containing donor IDs (Option B above). Default \code{NULL}.
#' @param donors Character vector or \code{NULL}. Specific donor IDs to
#'   include. If \code{NULL} (default), all donors are used.
#' @param min_donors Integer. Minimum number of donors an interaction must
#'   appear in. Default \code{2}.
#' @param p Numeric. RBP patience parameter (0 < p < 1). Default \code{0.8}.
#' @param resource Character. LIANA ligand-receptor resource. Default
#'   \code{"CellPhoneDB"}.
#' @param nperms Integer. Permutations for CellPhoneDB. Use \code{100} for
#'   quick testing, \code{1000} for publication. Default \code{1000}.
#' @param expr_prop Numeric. Minimum expression proportion for LIANA.
#'   Default \code{0.05}.
#' @param workers Integer. Parallel workers for CellPhoneDB. Default \code{4}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A data frame sorted by \code{total_rbp_score} (descending) with
#'   columns: \code{source}, \code{target}, \code{ligand}, \code{receptor},
#'   \code{LR}, \code{total_rbp_score}, \code{n_donors}, and per-method
#'   RBP scores.
#'
#' @examples
#' \dontrun{
#' library(CCIcompass)
#'
#' # Option A: donor ID embedded in cell-type label
#' obj <- readRDS("my_seurat_object.rds")
#' results <- run_pipeline(
#'   seurat_obj = obj,
#'   idents_col = "celltype.combo",  # values like "Epi-Intestine_donor1"
#'   min_donors = 2
#' )
#'
#' # Option B: separate donor column
#' results <- run_pipeline(
#'   seurat_obj = obj,
#'   idents_col = "celltype",        # values like "Epi-Intestine"
#'   donor_col  = "patient_id",      # values like "donor1"
#'   min_donors = 2
#' )
#'
#' plot_cis_dotplot(results, top_n = 20)
#' }
#'
#' @importFrom dplyr mutate filter group_by ungroup arrange desc
#'   n_distinct slice inner_join rename_with all_of
#' @export
run_pipeline <- function(seurat_obj,
                         idents_col = "celltype.combo",
                         donor_col  = NULL,
                         donors     = NULL,
                         min_donors = 2,
                         p          = 0.8,
                         resource   = "CellPhoneDB",
                         nperms     = 1000,
                         expr_prop  = 0.05,
                         workers    = 4,
                         verbose    = TRUE) {

  .msg <- function(...) if (verbose) message("[CCIcompass] ", ...)

  # ── Validate inputs ──────────────────────────────────────────────────────
  if (!inherits(seurat_obj, "Seurat"))
    stop("`seurat_obj` must be a Seurat object.", call. = FALSE)

  meta <- seurat_obj@meta.data

  if (!idents_col %in% colnames(meta))
    stop("`idents_col` '", idents_col, "' not found in meta.data.\n",
         "Available columns: ", paste(colnames(meta), collapse = ", "),
         call. = FALSE)

  if (!is.null(donor_col) && !donor_col %in% colnames(meta))
    stop("`donor_col` '", donor_col, "' not found in meta.data.", call. = FALSE)

  if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1)
    stop("`p` must be a single number strictly between 0 and 1.", call. = FALSE)

  if (!requireNamespace("liana", quietly = TRUE))
    stop("Package 'liana' is required. Install with:\n",
         "  devtools::install_github('saezlab/liana')", call. = FALSE)

  # ── Step 0: Build combined cell-type + donor label if needed ─────────────
  if (!is.null(donor_col)) {
    .msg("Building combined cell-type + donor label...")
    seurat_obj@meta.data[["ccicompass_idents"]] <-
      paste(meta[[idents_col]], meta[[donor_col]], sep = "_")
    run_col <- "ccicompass_idents"
  } else {
    run_col <- idents_col
  }

  # Auto-detect donors from labels
  if (is.null(donors)) {
    all_labels <- unique(seurat_obj@meta.data[[run_col]])
    donors <- unique(sub(".*_([^_]+)$", "\\1", all_labels))
    .msg("Detected ", length(donors), " donor(s): ",
         paste(donors, collapse = ", "))
  }

  # ── Step 1: Run LIANA (all 6 methods) ────────────────────────────────────
  .msg("Step 1/5: Running LIANA — 4 rank-based methods...")

  liana_multi <- liana::liana_wrap(
    seurat_obj,
    idents_col = run_col,
    method     = c("connectome", "logfc", "natmi", "sca"),
    resource   = resource
  )

  .msg("Step 1/5: Running CellPhoneDB (", nperms, " permutations)...")
  options(future.globals.maxSize = 10000 * 1024^2)
  liana_cpdb <- liana::liana_wrap(
    seurat_obj,
    idents_col         = run_col,
    method             = "cellphonedb",
    resource           = resource,
    permutation.params = list(nperms      = nperms,
                              parallelize = FALSE,
                              workers     = workers),
    expr_prop          = expr_prop
  )

  .msg("Step 1/5: Running CellChat...")
  liana_cc <- liana::liana_wrap(
    seurat_obj,
    idents_col = run_col,
    method     = "call_cellchat",
    resource   = resource
  )

  # Assemble named list and standardise column names
  method_list <- list(
    connectome  = liana_multi[["connectome"]],
    logfc       = liana_multi[["logfc"]],
    natmi       = liana_multi[["natmi"]],
    sca         = liana_multi[["sca"]],
    cellphonedb = as.data.frame(liana_cpdb),
    cellchat    = as.data.frame(liana_cc)
  )

  # Rename score columns to consistent names
  score_renames <- list(
    connectome  = c("weight_sc"  = "connectome_weight_sc"),
    logfc       = c("logfc_comb" = "logfc_logfc_comb"),
    natmi       = c("prod_weight"= "natmi_prod_weight"),
    sca         = c("LRscore"    = "sca_LRscore"),
    cellphonedb = c("lr.mean"    = "cellphonedb_lr.mean"),
    cellchat    = c("prob"       = "cellchat_prob")
  )

  method_list <- lapply(names(method_list), function(m) {
    df <- method_list[[m]]
    # Rename score column
    for (old in names(score_renames[[m]])) {
      new <- score_renames[[m]][[old]]
      if (old %in% colnames(df)) colnames(df)[colnames(df) == old] <- new
    }
    # liana uses ligand.complex / receptor.complex -> alias
    if ("ligand.complex"   %in% colnames(df)) df$ligand   <- df$ligand.complex
    if ("receptor.complex" %in% colnames(df)) df$receptor <- df$receptor.complex
    df
  })
  names(method_list) <- c("connectome","logfc","natmi","sca",
                           "cellphonedb","cellchat")

  # ── Step 2: Filter to same-donor pairs ───────────────────────────────────
  .msg("Step 2/5: Filtering to same-donor source-target pairs...")

  extract_donor <- function(x) sub(".*_([^_]+)$", "\\1", x)
  strip_donor   <- function(x) sub("_[^_]+$", "", x)

  filter_same_donor <- function(df) {
    df |>
      dplyr::mutate(
        .src_d     = extract_donor(source),
        .tgt_d     = extract_donor(target),
        patient_id = .src_d
      ) |>
      dplyr::filter(.src_d == .tgt_d, patient_id %in% donors) |>
      dplyr::mutate(source = strip_donor(source),
                    target = strip_donor(target)) |>
      dplyr::select(-.src_d, -.tgt_d)
  }

  method_list <- lapply(method_list, filter_same_donor)

  # ── Step 3: Keep interactions present in >= min_donors donors ────────────
  .msg("Step 3/5: Retaining interactions in >= ", min_donors, " donors...")

  keep_multi_donor <- function(df) {
    df |>
      dplyr::group_by(source, target, ligand, receptor) |>
      dplyr::filter(dplyr::n_distinct(patient_id) >= min_donors) |>
      dplyr::ungroup()
  }

  method_list <- lapply(method_list, keep_multi_donor)

  # ── Step 4: Intersect across all six methods ──────────────────────────────
  .msg("Step 4/5: Intersecting across all 6 methods...")

  join_keys <- c("patient_id", "source", "target", "ligand", "receptor")

  add_prefix <- function(df, prefix) {
    non_key <- setdiff(colnames(df), join_keys)
    dplyr::rename_with(df, ~ paste0(prefix, "_", .), dplyr::all_of(non_key))
  }

  prefixed   <- mapply(add_prefix, df = method_list,
                       prefix = names(method_list), SIMPLIFY = FALSE)
  consensus  <- prefixed[[1]]
  for (i in seq(2, length(prefixed)))
    consensus <- dplyr::inner_join(consensus, prefixed[[i]], by = join_keys)

  if (nrow(consensus) == 0L)
    stop("No interactions survived the 6-method consensus filter.\n",
         "Tips: lower `min_donors`, raise `expr_prop`, or check that all ",
         "6 methods produced non-empty output.", call. = FALSE)

  .msg("  Consensus interactions: ", nrow(consensus))

  # ── Step 5: Compute CIS ───────────────────────────────────────────────────
  .msg("Step 5/5: Computing Composite Interaction Score...")

  final <- compute_cis(consensus, p = p) |>
    dplyr::group_by(source, target, ligand, receptor) |>
    dplyr::mutate(n_donors = dplyr::n_distinct(patient_id)) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::mutate(LR = paste(ligand, receptor, sep = "-")) |>
    dplyr::arrange(dplyr::desc(total_rbp_score))

  .msg("Done. Returning ", nrow(final), " ranked interactions.")
  final
}
