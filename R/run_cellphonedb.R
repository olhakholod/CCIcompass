#' Run LIANA with CellPhoneDB-style permutation
#'
#' @param seurat_obj A Seurat object
#' @param idents_col Column for cell types
#' @return Data frame of inferred CCIs
#' @export
run_cellphonedb <- function(seurat_obj, idents_col = "celltype.combo") {
  liana::liana_wrap(
    seurat_obj,
    idents_col = idents_col,
    method = 'cellphonedb',
    resource = c('CellPhoneDB'),
    permutation.params = list(nperms = 100, parallelize = FALSE, workers = 4),
    expr_prop = 0.05
  )
}