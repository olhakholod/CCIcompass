#' Calculate Composite Interaction Score (CIS)
#'
#' @param df Data frame with common interactions
#' @param p Patience parameter (default = 0.9)
#' @return Data frame with CIS scores and rankings
#' @export
calculate_CIS <- function(df, p = 0.9) {
  methods_columns <- list(
    cellchat = "cellchat_prob",
    cellphonedb = "cellphonedb_lr.mean",
    connectome = "connectome_weight_sc",
    logfc = "logfc_logfc_comb",
    natmi = "natmi_prod_weight",
    sca = "sca_LRscore"
  )

  # Rank and compute RBP
  for (method in names(methods_columns)) {
    col <- methods_columns[[method]]
    rank_col <- paste0(method, "_rank")
    rbp_col <- paste0("rbp_", method)
    df[[rank_col]] <- rank(-df[[col]], ties.method = "min")
    df[[rbp_col]] <- p^(df[[rank_col]] - 1)
  }

  # Sum RBP scores
  df$CIS <- rowSums(df %>% select(starts_with("rbp_")))

  df <- df %>% arrange(desc(CIS))
  return(df)
}
