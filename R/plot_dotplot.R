#' Dot plot of top-ranked CCIs by CIS
#'
#' Creates a ggplot2 dot plot showing ligand-receptor pairs on the y-axis
#' and source-target cell pairs on the x-axis. Point size encodes CIS and
#' point colour encodes tissue of origin (if a \code{tissue} column exists).
#'
#' @param data A data frame from \code{\link{run_pipeline}} or
#'   \code{\link{compute_cis}}. Must contain \code{source}, \code{target},
#'   \code{ligand}, \code{receptor}, and \code{total_rbp_score}.
#' @param top_n Integer. Number of top LR pairs to show (ranked by maximum
#'   CIS across all cell pairs). Default \code{20}.
#' @param tissue_colours Named character vector mapping tissue names to
#'   hex colours. Only used when \code{data} contains a \code{tissue} column.
#'   Default uses orange/blue/green for intestine/skin/uterus.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param point_size_range Numeric vector of length 2. Min and max point
#'   sizes. Default \code{c(1, 7)}.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data("example_consensus")
#' cis <- compute_cis(example_consensus)
#' plot_cis_dotplot(cis, top_n = 10, title = "Top 10 CCIs")
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous
#'   scale_colour_manual coord_flip theme_bw labs theme element_text
#'   element_line unit
#' @importFrom dplyr mutate filter group_by summarise arrange desc
#'   slice_head pull
#' @export
plot_cis_dotplot <- function(data,
                              top_n            = 20,
                              tissue_colours   = c(intestine = "#E5835A",
                                                   skin      = "#82B0D9",
                                                   uterus    = "#7DC99E"),
                              title            = NULL,
                              point_size_range = c(1, 7)) {

  if (!is.data.frame(data))
    stop("`data` must be a data frame.", call. = FALSE)

  required <- c("source", "target", "ligand", "receptor", "total_rbp_score")
  missing  <- setdiff(required, colnames(data))
  if (length(missing) > 0L)
    stop("Missing columns: ", paste(missing, collapse = ", "), call. = FALSE)

  data <- data |>
    dplyr::mutate(
      LR            = paste(ligand, receptor, sep = "-"),
      source_target = paste(source, target,   sep = "-")
    )

  top_lrs <- data |>
    dplyr::group_by(LR) |>
    dplyr::summarise(max_cis = max(total_rbp_score, na.rm = TRUE),
                     .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_cis)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull(LR)

  plot_data <- data |>
    dplyr::filter(LR %in% top_lrs) |>
    dplyr::mutate(LR = factor(LR, levels = rev(sort(top_lrs))))

  p <- ggplot2::ggplot(plot_data,
         ggplot2::aes(x = LR, y = source_target,
                      size = total_rbp_score)) +
    ggplot2::scale_size_continuous(
      name   = "CIS",
      range  = point_size_range,
      breaks = seq(0.5, 5, 0.5)
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1,
                                               vjust = 1, size = 7),
      axis.text.y     = ggplot2::element_text(size = 7),
      panel.grid      = ggplot2::element_line(colour = "grey92"),
      legend.key.size = ggplot2::unit(0.4, "cm")
    ) +
    ggplot2::labs(x = "CCIs", y = "source-target cells", title = title)

  if ("tissue" %in% colnames(plot_data)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = tissue), alpha = 0.85) +
      ggplot2::scale_colour_manual(values = tissue_colours, name = NULL)
  } else {
    p <- p + ggplot2::geom_point(colour = "#4A90A4", alpha = 0.85)
  }

  p
}
