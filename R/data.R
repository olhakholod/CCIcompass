#' Example CCI consensus table
#'
#' A data frame of cell-cell interactions retained after intersecting results
#' across all six CCI inference methods and filtering to interactions common
#' across all three donors. This is the direct input to
#' \code{\link{compute_cis}}.
#'
#' @format A data frame with 640 rows and 11 columns:
#' \describe{
#'   \item{patient_id}{Donor identifier (character).}
#'   \item{source}{Sending cell type (character).}
#'   \item{target}{Receiving cell type (character).}
#'   \item{ligand}{Ligand gene symbol (character).}
#'   \item{receptor}{Receptor gene symbol (character).}
#'   \item{cellchat_prob}{CellChat interaction probability (numeric).}
#'   \item{cellphonedb_lr.mean}{CellPhoneDB mean LR expression (numeric).}
#'   \item{connectome_weight_sc}{Connectome edge weight (numeric).}
#'   \item{logfc_logfc_comb}{LogFC combined score (numeric).}
#'   \item{natmi_prod_weight}{NATMI product weight (numeric).}
#'   \item{sca_LRscore}{SCA ligand-receptor score (numeric).}
#' }
#' @source Synthetic data matching the structure of LIANA outputs from:
#'   Kholod O. et al., PLOS Computational Biology (2025).
"example_consensus"


#' Example per-method CCI list
#'
#' A named list of six data frames, each containing the output of one CCI
#' inference method as produced by \code{liana_wrap()} on a small subset
#' of intestinal single-cell RNA-seq data.
#'
#' @format A named list with elements \code{cellchat}, \code{cellphonedb},
#'   \code{connectome}, \code{logfc}, \code{natmi}, and \code{sca}.
#'   Each data frame has columns \code{patient_id}, \code{source},
#'   \code{target}, \code{ligand}, \code{receptor}, and one method-specific
#'   score column.
#' @source Synthetic data. See \code{\link{example_consensus}}.
"example_cci_list"


#' HPA MIF immunohistochemistry data
#'
#' Percentage of MIF-positive cells per donor and cell type, scored from
#' Human Protein Atlas (HPA) immunohistochemistry images.
#'
#' @format A data frame with 9 rows and 4 columns:
#' \describe{
#'   \item{donor}{Donor label, e.g. \code{"donor 1"} (character).}
#'   \item{cell_type}{Cell type: one of \code{"Glandular cells of intestine"},
#'     \code{"Keratinocytes of skin"}, or \code{"Glandular cells of uterus"}.}
#'   \item{pct_positive}{Mean percentage of MIF-positive cells (numeric).}
#'   \item{sd}{Standard deviation across regions of interest (numeric).}
#' }
#' @source Kholod O. et al., PLOS Computational Biology (2025).
"example_hpa_mif"


#' HPA CD74 immunohistochemistry data
#'
#' Percentage of CD74-positive cells per donor and cell type, scored from
#' Human Protein Atlas (HPA) immunohistochemistry images.
#'
#' @format A data frame with 6 rows and 4 columns:
#' \describe{
#'   \item{donor}{Donor label (character).}
#'   \item{cell_type}{Cell type: one of \code{"Glandular cells of intestine"}
#'     or \code{"Glandular cells of uterus"}.}
#'   \item{pct_positive}{Mean percentage of CD74-positive cells (numeric).}
#'   \item{sd}{Standard deviation across regions of interest (numeric).}
#' }
#' @source Kholod O. et al., PLOS Computational Biology (2025).
"example_hpa_cd74"
