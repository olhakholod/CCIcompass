#' Run LIANA with multiple methods
#'
#' @param seurat_obj A Seurat object
#' @param idents_col Column name in metadata for identities
#' @return A named list of data frames from different methods
#' @export
run_liana <- function(seurat_obj, idents_col = "celltype.combo") {
  liana::liana_wrap(
    seurat_obj,
    idents_col = idents_col,
    method = c("connectome", "logfc", "natmi", "sca"),
    resource = c("CellPhoneDB")
  )
}