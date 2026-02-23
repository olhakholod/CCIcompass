#' Compute the Composite Interaction Score (CIS)
#'
#' Ranks cell–cell interactions from six CCI inference methods using
#' Rank-Biased Precision (RBP) and sums the scores into a single
#' Composite Interaction Score (CIS) per interaction.
#'
#' @details
#' For each method \eqn{m}, an interaction at rank \eqn{r} receives an RBP
#' score of \eqn{p^{r-1}}, where \eqn{p} is the patience parameter
#' (0 < p < 1). A value of \eqn{p = 0.8} means the top-ranked interaction
#' gets a score of 1.0, the second gets 0.8, the third 0.64, and so on.
#' The CIS is the sum of RBP scores across all six methods, giving a maximum
#' possible value of 6.0 (interaction ranked first by every method).
#'
#' Expected input columns (produced by \code{\link{intersect_methods}}):
#' \itemize{
#'   \item \code{patient_id}, \code{source}, \code{target},
#'         \code{ligand}, \code{receptor}
#'   \item \code{cellchat_prob}
#'   \item \code{cellphonedb_lr.mean}
#'   \item \code{connectome_weight_sc}
#'   \item \code{logfc_logfc_comb}
#'   \item \code{natmi_prod_weight}
#'   \item \code{sca_LRscore}
#' }
#'
#' @param data A data frame of consensus CCIs with per-method score columns
#'   (output of \code{\link{intersect_methods}}).
#' @param p Numeric scalar. Patience parameter for RBP (default \code{0.8}).
#'   Must satisfy \code{0 < p < 1}. Higher values weight top-ranked
#'   interactions more strongly.
#'
#' @return The input data frame with the following additional columns:
#' \describe{
#'   \item{\code{*_rank}}{Integer rank per method (1 = highest score).}
#'   \item{\code{rbp_*}}{RBP score per method.}
#'   \item{\code{total_rbp_score}}{CIS: sum of RBP scores across all methods.}
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
  if (!is.numeric(p) || length(p) != 1 || p <= 0 || p >= 1)
    stop("`p` must be a single numeric value strictly between 0 and 1.",
         call. = FALSE)

  required_cols <- c(
    "cellchat_prob", "cellphonedb_lr.mean", "connectome_weight_sc",
    "logfc_logfc_comb", "natmi_prod_weight", "sca_LRscore"
  )
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0)
    stop("Missing required columns: ", paste(missing, collapse = ", "),
         call. = FALSE)

  rbp <- function(rank_vec, patience) patience ^ (rank_vec - 1)

  data <- data |>
    dplyr::mutate(
      cellchat_rank    = rank(-cellchat_prob,        ties.method = "min"),
      cellphonedb_rank = rank(-cellphonedb_lr.mean,  ties.method = "min"),
      connectome_rank  = rank(-connectome_weight_sc, ties.method = "min"),
      logfc_rank       = rank(-logfc_logfc_comb,     ties.method = "min"),
      natmi_rank       = rank(-natmi_prod_weight,    ties.method = "min"),
      sca_rank         = rank(-sca_LRscore,          ties.method = "min")
    ) |>
    dplyr::mutate(
      rbp_cellchat    = rbp(cellchat_rank,    p),
      rbp_cellphonedb = rbp(cellphonedb_rank, p),
      rbp_connectome  = rbp(connectome_rank,  p),
      rbp_logfc       = rbp(logfc_rank,       p),
      rbp_natmi       = rbp(natmi_rank,       p),
      rbp_sca         = rbp(sca_rank,         p)
    ) |>
    dplyr::mutate(
      total_rbp_score = rbp_cellchat + rbp_cellphonedb + rbp_connectome +
                        rbp_logfc    + rbp_natmi        + rbp_sca
    ) |>
    dplyr::arrange(dplyr::desc(total_rbp_score))

  data
}
