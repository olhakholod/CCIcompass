#' Dot plot of top-ranked CCIs by CIS
#'
#' Creates a ggplot2 dot plot showing ligand–receptor pairs (y-axis) against
#' source–target cell pairs (x-axis). Point size encodes the CIS value and
#' point colour encodes tissue of origin.
#'
#' @param data A data frame containing CCI results. Must include columns for
#'   \code{source}, \code{target}, \code{ligand}, \code{receptor}, and
#'   \code{total_rbp_score}. An optional \code{tissue} column provides the
#'   colour aesthetic.
#' @param top_n Integer. Number of top-ranked LR pairs to display (ranked by
#'   maximum \code{total_rbp_score} across source–target pairs).
#'   Default \code{20}.
#' @param tissue_colours Named character vector mapping tissue names to hex
#'   colours. Default uses orange/blue/green for intestine/skin/uterus.
#' @param title Character. Plot title. Default \code{NULL}.
#' @param x_angle Numeric. Angle of x-axis text labels. Default \code{45}.
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
#' @importFrom dplyr mutate filter group_by summarise arrange desc slice_head
#' @export
plot_cis_dotplot <- function(data,
                              top_n          = 20,
                              tissue_colours = c(intestine = "#E5835A",
                                                 skin      = "#82B0D9",
                                                 uterus    = "#7DC99E"),
                              title          = NULL,
                              x_angle        = 45) {

  data <- data |>
    dplyr::mutate(
      LR           = paste(ligand, receptor, sep = "-"),
      source_target = paste(source, target, sep = "-")
    )

  # Select top_n LR pairs by max CIS across all cell pairs
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
    ggplot2::geom_point(alpha = 0.85,
                        colour = if ("tissue" %in% names(plot_data))
                          NULL else "#E5835A") +
    ggplot2::scale_size_continuous(
      name   = "CIS",
      range  = c(1, 7),
      breaks = seq(0.5, 4.5, 0.5)
    ) +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(angle = x_angle, hjust = 1,
                                            vjust = 1, size = 7),
      axis.text.y  = ggplot2::element_text(size = 7),
      panel.grid   = ggplot2::element_line(colour = "grey92"),
      legend.key.size = ggplot2::unit(0.4, "cm")
    ) +
    ggplot2::labs(x = "CCIs", y = "source–target cells", title = title)

  if ("tissue" %in% names(plot_data)) {
    p <- p +
      ggplot2::aes(colour = tissue) +
      ggplot2::scale_colour_manual(values = tissue_colours, name = NULL)
  }

  p
}
