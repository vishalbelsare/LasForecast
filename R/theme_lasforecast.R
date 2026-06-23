#' Publication-ready ggplot2 theme for LasForecast
#'
#' A clean theme designed for academic papers and presentations.
#' Used consistently across all package visualization functions.
#'
#' @param base_size Base font size (default 11).
#' @param base_family Base font family (default "").
#'
#' @return A ggplot2 theme object.
#'
#' @import ggplot2
#' @export
theme_lasforecast <- function(base_size = 11, base_family = "") {
    theme(
        # Panel
        panel.background = element_blank(),
        panel.border = element_rect(linetype = 1, colour = "black", fill = NA),
        panel.grid.major = element_line(linetype = 2, color = "grey90"),
        panel.grid.minor = element_blank(),

        # Strip (facets)
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = base_size),

        # Axis
        axis.text = element_text(size = base_size * 0.9, color = "grey30"),
        axis.title = element_text(size = base_size),
        axis.ticks = element_line(color = "grey70", linewidth = 0.3),

        # Legend
        legend.position = "bottom",
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size = base_size * 0.9),

        # Title
        plot.title = element_text(size = base_size * 1.2, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = base_size * 0.95, color = "grey40", hjust = 0.5),

        # Margin
        plot.margin = margin(10, 10, 10, 10)
    )
}


#' Generate ggplot2 default color palette
#'
#' Replicates the default ggplot2 color palette (equally spaced hues).
#'
#' @param n Number of colors.
#' @return Character vector of hex color codes.
#' @keywords internal
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    grDevices::hcl(h = hues, l = 65, c = 100)[1:n]
}
