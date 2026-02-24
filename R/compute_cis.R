#' Compute the Composite Interaction Score (CIS)
#'
#' Ranks cell-cell interactions from six CCI inference methods using
#' Rank-Biased Precision (RBP) and sums the scores into a single
#' Composite Interaction Score (CIS) per interaction.
#'
#' @details
#' For each method \eqn{m}, an interaction at rank \eqn{r} receives an RBP
#' score of \eqn{p^{r-1}}, where \eqn{p} is the patience parameter
#' (0 < p < 1). The CIS is the sum of RBP scores across all six methods,
#' giving a maximum possible value of 6.0 (ranked first by every method).
#'
#' Required columns (produced by \code{\link{run_pipeline}} or by manually
#' joining LIANA outputs):
#' \itemize{
#'   \item \code{cellchat_prob}
#'   \item \code{cellphonedb_lr.mean}
#'   \item \code{connectome_weight_sc}
#'   \item \code{logfc_logfc_comb}
#'   \item \code{natmi_prod_weight}
#'   \item \code{sca_LRscore}
#' }
#'
#' @param data A data frame with the six per-method score columns listed above.
#' @param p Numeric. Patience parameter (0 < p < 1). Default \code{0.8}.
#'
#' @return The input data frame with added columns:
#' \describe{
#'   \item{\code{*_rank}}{Integer rank per method (1 = best score).}
#'   \item{\code{rbp_*}}{RBP score per method.}
#'   \item{\code{total_rbp_score}}{CIS: sum of all six RBP scores.}
#' }
#' Rows are sorted by \code{total_rbp_score} descending.
#'
#' @examples
#' data("example_consensus")
#' cis <- compute_cis(example_consensus, p = 0.8)
#' head(cis[, c("ligand", "receptor", "source", "target", "total_rbp_score")])
#'
#' @importFrom dplyr mutate arrange desc
#' @export
compute_cis <- function(data, p = 0.8) {

  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)
  if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1)
    stop("`p` must be a single number strictly between 0 and 1.", call. = FALSE)

  required <- c("cellchat_prob", "cellphonedb_lr.mean", "connectome_weight_sc",
                "logfc_logfc_comb", "natmi_prod_weight", "sca_LRscore")
  missing  <- setdiff(required, colnames(data))
  if (length(missing) > 0L)
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)

  rbp <- function(rank_vec, patience) patience ^ (rank_vec - 1L)

  data |>
    dplyr::mutate(
      cellchat_rank    = rank(-cellchat_prob,        ties.method = "min"),
      cellphonedb_rank = rank(-cellphonedb_lr.mean,  ties.method = "min"),
      connectome_rank  = rank(-connectome_weight_sc, ties.method = "min"),
      logfc_rank       = rank(-logfc_logfc_comb,     ties.method = "min"),
      natmi_rank       = rank(-natmi_prod_weight,    ties.method = "min"),
      sca_rank         = rank(-sca_LRscore,          ties.method = "min"),
      rbp_cellchat     = rbp(cellchat_rank,    p),
      rbp_cellphonedb  = rbp(cellphonedb_rank, p),
      rbp_connectome   = rbp(connectome_rank,  p),
      rbp_logfc        = rbp(logfc_rank,       p),
      rbp_natmi        = rbp(natmi_rank,       p),
      rbp_sca          = rbp(sca_rank,         p),
      total_rbp_score  = rbp_cellchat + rbp_cellphonedb + rbp_connectome +
                         rbp_logfc    + rbp_natmi        + rbp_sca
    ) |>
    dplyr::arrange(dplyr::desc(total_rbp_score))
}
