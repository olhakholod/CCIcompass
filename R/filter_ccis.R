#' Filter CCIs to same-patient source–target pairs
#'
#' Removes interactions where the source and target cell types come from
#' different donors. This is necessary when a Seurat object contains cells
#' from multiple donors whose IDs are embedded in the cell-type label
#' (e.g. \code{"Epi-Intestine_A26"}).
#'
#' The function extracts the donor ID from each label as the substring
#' following the last underscore, then keeps only rows where the
#' source and target donor IDs match.
#'
#' @param data A data frame with at least \code{source} and \code{target}
#'   columns whose values encode the donor ID after an underscore
#'   (e.g. \code{"Epi-Intestine_A26 (386C)"}).
#' @param out_path Optional character. If provided, the filtered data frame
#'   is written to this CSV path. Default \code{NULL} (no file written).
#'
#' @return A data frame with only same-donor interactions.
#'
#' @examples
#' df <- data.frame(
#'   source = c("Epi_donor1", "Epi_donor1", "Epi_donor2"),
#'   target = c("Mac_donor1", "Mac_donor2", "Mac_donor2"),
#'   ligand = c("MIF",  "MIF",  "APP"),
#'   receptor = c("CD74", "CD74", "CD74")
#' )
#' filter_same_patient(df)
#'
#' @importFrom dplyr mutate filter select
#' @export
filter_same_patient <- function(data, out_path = NULL) {

  extract_id <- function(x) sub(".*_(.*?)\\s*$", "\\1", x)

  filtered <- data |>
    dplyr::mutate(
      .src_id = extract_id(source),
      .tgt_id = extract_id(target)
    ) |>
    dplyr::filter(.src_id == .tgt_id) |>
    dplyr::select(-.src_id, -.tgt_id)

  if (!is.null(out_path)) {
    utils::write.csv(filtered, out_path, row.names = FALSE)
    message("Filtered data written to: ", out_path)
  }

  filtered
}


#' Retain CCIs common across all (or a specified number of) donors
#'
#' Filters a CCI data frame to keep only ligand–receptor + source–target
#' combinations that appear in the required number of donors.
#'
#' @param data A data frame with columns \code{source}, \code{target},
#'   \code{ligand}, \code{receptor}, and \code{patient_id}.
#' @param donors Character vector of donor IDs to include.
#' @param n_req Integer. Minimum number of donors in which an interaction must
#'   appear. Defaults to \code{length(donors)} (i.e. all donors).
#' @param out_path Optional CSV output path. Default \code{NULL}.
#'
#' @return A filtered data frame.
#'
#' @examples
#' data("example_cci_list")
#' df <- example_cci_list$sca
#' filtered <- filter_common_donors(
#'   df,
#'   donors = unique(df$patient_id),
#'   n_req  = 2
#' )
#'
#' @importFrom dplyr filter group_by ungroup
#' @export
filter_common_donors <- function(data, donors, n_req = length(donors),
                                  out_path = NULL) {

  filtered <- data |>
    dplyr::filter(patient_id %in% donors) |>
    dplyr::group_by(source, target, ligand, receptor) |>
    dplyr::filter(dplyr::n() >= n_req) |>
    dplyr::ungroup()

  if (!is.null(out_path)) {
    utils::write.csv(filtered, out_path, row.names = FALSE)
    message("Common-donor filtered data written to: ", out_path)
  }

  filtered
}


#' Intersect CCIs across multiple inference methods
#'
#' Performs an inner join across a named list of per-method data frames,
#' retaining only interactions present in **all** methods.
#' Score columns are renamed with a method prefix to avoid conflicts.
#'
#' @param method_list A named list of data frames, one per method. Each must
#'   contain \code{patient_id}, \code{source}, \code{target}, \code{ligand},
#'   \code{receptor} as key columns. Additional columns are the method scores.
#' @param out_path Optional CSV output path. Default \code{NULL}.
#'
#' @return A single data frame containing only consensus interactions, with
#'   all method score columns present (prefixed by method name).
#'
#' @examples
#' \dontrun{
#' consensus <- intersect_methods(
#'   list(sca = sca_df, natmi = natmi_df, logfc = logfc_df)
#' )
#' }
#'
#' @importFrom dplyr rename_with inner_join
#' @export
intersect_methods <- function(method_list, out_path = NULL) {

  if (length(method_list) < 2)
    stop("`method_list` must contain at least two methods.", call. = FALSE)

  key_cols <- c("patient_id", "source", "target", "ligand", "receptor")

  prefixed <- lapply(names(method_list), function(m) {
    df <- method_list[[m]]
    dplyr::rename_with(df, ~paste0(m, "_", .), -dplyr::any_of(key_cols))
  })

  common <- prefixed[[1]]
  for (i in seq(2, length(prefixed))) {
    common <- dplyr::inner_join(common, prefixed[[i]], by = key_cols)
  }

  if (!is.null(out_path)) {
    utils::write.csv(common, out_path, row.names = FALSE)
    message("Consensus interactions written to: ", out_path)
  }

  common
}


#' Filter interactions to those present in a minimum number of donors
#'
#' After CIS computation, applies a final quality filter to keep only
#' interactions observed in at least \code{min_donors} independent donors.
#' One representative row per interaction is retained (the first occurrence,
#' which has the highest CIS after \code{\link{compute_cis}} sorting).
#'
#' @param data A data frame with columns \code{source}, \code{target},
#'   \code{ligand}, \code{receptor}, and \code{patient_id}.
#' @param min_donors Integer. Minimum number of distinct donors required.
#'   Default \code{2}.
#' @param out_path Optional CSV output path. Default \code{NULL}.
#'
#' @return A filtered data frame with one row per unique interaction.
#'
#' @examples
#' data("example_consensus")
#' cis <- compute_cis(example_consensus)
#' final <- filter_min_donors(cis, min_donors = 2)
#'
#' @importFrom dplyr group_by filter slice ungroup n_distinct
#' @export
filter_min_donors <- function(data, min_donors = 2, out_path = NULL) {

  filtered <- data |>
    dplyr::group_by(source, target, ligand, receptor) |>
    dplyr::filter(dplyr::n_distinct(patient_id) >= min_donors) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  if (!is.null(out_path)) {
    utils::write.csv(filtered, out_path, row.names = FALSE)
    message("Final filtered interactions written to: ", out_path)
  }

  filtered
}
