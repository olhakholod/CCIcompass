#' Run LIANA with CellChat method
#'
#' @param seurat_obj A Seurat object
#' @param idents_col Column for cell types
#' @return Data frame of inferred CCIs
#' @export
run_cellchat <- function(seurat_obj, idents_col = "celltype.combo") {
  liana::liana_wrap(
    seurat_obj,
    idents_col = idents_col,
    method = c("call_cellchat"),
    resource = c("CellPhoneDB")
  )
}