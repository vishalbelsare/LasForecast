#' Trend plot
#' @param y_0 True predict target, length n vector
#' @param y_hat Predicted values by m different methods, n-by-m matrix, colnames well defined
#' @param dates Date vector of class "Date". \cr
#' Use lubridate and zoo to transfer original date vector to "Date" class
#' @param pt_num Number date points (pt_num + 1) wanted on the x-axis
#' @param col_vec A vector indicates colors of different trends. \cr
#' e.g.:  c("#D8DBE2", "#F46B7B", "#518DE8", "#FFBC42" )\cr
#' If NULL, use ggplot default.
#' @param line_size line size for each trend
#' @param alpha_size degree of appearance
#' @param xlab x-axis label
#' @param ylab y-axis label
#'
#' @return A ggplot output
#'
#' @export plot_trend
plot_trend <- function(y_0,
                       y_hat,
                       dates,
                       pt_num = 4,
                       col_vec = NULL,
                       line_size = rep(0.75, ncol(y_hat) + 1),
                       alpha_size = rep(0.75, ncol(y_hat) + 1),
                       xlab = NULL,
                       ylab = NULL) {
    n <- length(y_0)
    m <- ncol(y_hat)
    methods_use <- colnames(y_hat)

    if (!inherits(dates, "Date")) {
        dates <- zoo::as.Date(dates)
    }

    if (is.null(col_vec)) {
        col_vec <- gg_color_hue(m + 1)
    }

    # Prepare data frame
    df_plot <- cbind(data.frame(date = dates, y_0 = y_0), as.data.frame(y_hat))
    df_plot <- tidyr::pivot_longer(
        df_plot,
        cols = -"date",
        names_to = "variable",
        values_to = "value"
    )

    pt <- dates[seq(from = 1, to = n, by = max(1, floor(n / pt_num - 1)))]

    p_out <- ggplot(data = df_plot) +
        geom_line(mapping = aes(
            x = .data[["date"]],
            y = .data[["value"]],
            color = .data[["variable"]],
            linewidth = .data[["variable"]],
            alpha = .data[["variable"]]
        )) +
        scale_colour_manual(
            breaks = c("y_0", methods_use),
            values = col_vec,
            labels = c("True value", methods_use),
            guide = guide_legend(override.aes = aes(alpha = NA))
        ) +
        scale_linewidth_manual(
            breaks = c("y_0", methods_use),
            values = line_size,
            labels = c("True value", methods_use)
        ) +
        scale_alpha_manual(
            breaks = c("y_0", methods_use),
            values = alpha_size,
            labels = c("True value", methods_use)
        ) +
        scale_x_date(breaks = pt, labels = as.character(lubridate::year(pt))) +
        labs(x = xlab, y = ylab) +
        theme_lasforecast()

    return(p_out)
}

#' Coefficient plot
#'
#' @param coef_est estimated slope by m different methods, list of length m with each element n-by-p matrix\cr
#' names of list well defined
#' @param dates Date vector of class "Date".
#' @param pt_num Number date points (pt_num + 1) wanted on the x-axis
#' @param col_vec A vector indicates colors of different trends.
#' @param line_size line size for each trend
#' @param alpha_size degree of appearance
#' @param xlab x-axis label
#' @param ylab y-axis label
#' @param num_col number of columns in facet wrap
#'
#' @return A ggplot output
#'
#' @export plot_coef
plot_coef <- function(coef_est, dates,
                      pt_num = 4,
                      col_vec = NULL,
                      line_size = rep(0.75, length(coef_est)),
                      alpha_size = rep(0.75, length(coef_est)),
                      xlab = NULL,
                      ylab = NULL,
                      num_col = 1) {
    n <- nrow(coef_est[[1]])
    p <- ncol(coef_est[[1]])
    m <- length(coef_est)

    methods_use <- names(coef_est)
    var_names <- colnames(coef_est[[1]])

    if (!inherits(dates, "Date")) {
        dates <- zoo::as.Date(dates)
    }

    if (is.null(col_vec)) {
        col_vec <- gg_color_hue(m)
    }

    pt <- dates[seq(from = 1, to = n, by = max(1, floor(n / pt_num - 1)))]

    df_plot <- NULL
    for (i in seq_len(m)) {
        method_i <- methods_use[i]
        df_temp <- cbind(data.frame(date = dates), as.data.frame(coef_est[[method_i]]))
        df_temp$method <- method_i
        df_plot <- rbind(df_plot, df_temp)
    }

    df_plot <- tidyr::pivot_longer(
        df_plot,
        cols = -c("date", "method"),
        names_to = "variable",
        values_to = "value"
    )
    df_plot$variable <- factor(df_plot$variable, levels = var_names)

    p_out <- ggplot(data = df_plot) +
        geom_line(mapping = aes(
            x = .data[["date"]],
            y = .data[["value"]],
            color = .data[["method"]],
            linewidth = .data[["method"]],
            alpha = .data[["method"]]
        )) +
        scale_x_date(breaks = pt, labels = as.character(lubridate::year(pt))) +
        labs(x = xlab, y = ylab) +
        theme_lasforecast() +
        scale_color_manual(
            values = col_vec,
            guide = guide_legend(override.aes = aes(alpha = NA))
        ) +
        scale_linewidth_manual(values = line_size) +
        scale_alpha_manual(values = alpha_size) +
        facet_wrap(~variable, ncol = num_col, scales = "free", strip.position = "right")

    return(p_out)
}

# Re-export the ggplot2 autoplot generic so autoplot() resolves to the methods
# below when only LasForecast is attached (without library(ggplot2)).
#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot


#' Autoplot method for lasforecast_roll objects
#'
#' @param object An object of class "lasforecast_roll" generated by \code{\link{roll_predict}}.
#' @param methods_to_plot Character vector of methods to plot. Default uses all.
#' @param date_range Length-2 vector of dates (or strings coercible to Date)
#'   to restrict the plotted time span. Ignored when dates are unavailable.
#' @param ... Additional arguments.
#'
#' @import ggplot2
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter select mutate
#' @importFrom stats setNames
#' @export
autoplot.lasforecast_roll <- function(object, methods_to_plot = NULL,
                                      date_range = NULL, ...) {
    y <- object$y
    h <- object$h
    dates <- object$dates

    methods_avail <- object$methods_use
    if (is.null(methods_to_plot)) {
        methods_to_plot <- methods_avail
    } else {
        methods_to_plot <- intersect(methods_to_plot, methods_avail)
    }

    num_forecast <- length(object[[methods_avail[1]]]$y_hat)
    has_dates <- !is.null(dates) && length(dates) == num_forecast

    df_plot <- data.frame(
        Index = seq_len(num_forecast),
        Actual = tail(y, num_forecast)
    )
    if (has_dates) df_plot$Date <- dates

    for (m in methods_to_plot) {
        df_plot[[m]] <- object[[m]]$y_hat
    }

    x_var <- if (has_dates) "Date" else "Index"

    if (!is.null(date_range) && has_dates) {
        date_range <- as.Date(date_range)
        df_plot <- df_plot[df_plot$Date >= date_range[1] &
                           df_plot$Date <= date_range[2], , drop = FALSE]
    }

    series <- c("Actual", methods_to_plot)
    df_long <- tidyr::pivot_longer(
        df_plot, cols = dplyr::all_of(series),
        names_to = "Series", values_to = "Value"
    )
    df_long$Series <- factor(df_long$Series, levels = series)

    method_colors <- gg_color_hue(length(methods_to_plot))
    colors <- c("Actual" = "grey50", setNames(method_colors, methods_to_plot))
    widths <- c("Actual" = 0.8, setNames(rep(0.5, length(methods_to_plot)), methods_to_plot))
    alphas <- c("Actual" = 0.4, setNames(rep(0.9, length(methods_to_plot)), methods_to_plot))

    p <- ggplot(df_long, aes(x = .data[[x_var]], y = .data[["Value"]],
                             color = .data[["Series"]],
                             linewidth = .data[["Series"]],
                             alpha = .data[["Series"]])) +
        geom_line() +
        scale_color_manual(values = colors) +
        scale_linewidth_manual(values = widths, guide = "none") +
        scale_alpha_manual(values = alphas, guide = "none") +
        theme_lasforecast() +
        labs(
            x = if (has_dates) "Date" else "Time Index",
            y = "Value",
            color = NULL
        )
    return(p)
}


#' Autoplot method for lasforecast_backtest objects
#'
#' @param object An object of class \code{"lasforecast_backtest"}.
#' @param type Type of plot. Loss-metric plots: \code{"bar"}, \code{"heatmap"}.
#'   Forecast plot: \code{"forecasts"}. Variable selection plots:
#'   \code{"frequency"} (selection frequency bars), \code{"coef_heatmap"}
#'   (time x variable coefficient heatmap), \code{"coef_path"} (coefficient
#'   evolution over time).
#' @param methods Character vector of methods to include. Default uses all.
#' @param loss Character vector of loss metrics to include (e.g., \code{"RMSE"},
#'   \code{"MAE"}). Used by \code{"bar"} and \code{"heatmap"} types.
#' @param ratio Logical. If \code{TRUE}, plot ratios relative to the benchmark
#'   instead of raw metric values. Default \code{FALSE}.
#' @param top_n Number of top variables to show for variable selection plots
#'   (default 20).
#' @param standardize Logical, forwarded to \code{\link{plot_variable_importance}}
#'   for the variable-selection plot types (\code{"frequency"},
#'   \code{"coef_heatmap"}, \code{"coef_path"}): if \code{TRUE}, each coefficient
#'   is scaled by its predictor's standard deviation. Default \code{NULL} uses
#'   the per-type default.
#' @param ... Additional arguments passed through (e.g., \code{methods_to_plot}
#'   and \code{date_range} for \code{type = "forecasts"}).
#'
#' @import ggplot2
#' @export
autoplot.lasforecast_backtest <- function(object,
                                          type = c("bar", "heatmap", "forecasts",
                                                   "frequency", "coef_heatmap",
                                                   "coef_path"),
                                          methods = NULL, loss = NULL,
                                          ratio = FALSE, top_n = 20,
                                          standardize = NULL, ...) {
    type <- match.arg(type)

    if (type == "forecasts") {
        return(autoplot(object$results, ...))
    }

    if (type %in% c("frequency", "coef_heatmap", "coef_path")) {
        vi_type <- switch(type,
            frequency = "frequency",
            coef_heatmap = "heatmap",
            coef_path = "evolution"
        )
        dots <- list(...)
        return(plot_variable_importance(object$results, methods = methods,
                                        top_n = top_n, type = vi_type,
                                        standardize = standardize,
                                        date_range = dots$date_range,
                                        variables = dots$variables))
    }

    # Build summary table with ratios
    st <- summary(object)
    benchmark <- object$benchmark

    # Identify raw loss columns (not _Ratio)
    raw_cols <- setdiff(names(st), c("Method", grep("_Ratio$", names(st), value = TRUE)))

    # Filter methods
    if (!is.null(methods)) {
        st <- st[st$Method %in% methods, , drop = FALSE]
    }

    # Select metric columns
    if (ratio) {
        metric_cols <- paste0(raw_cols, "_Ratio")
        metric_cols <- intersect(metric_cols, names(st))
    } else {
        metric_cols <- raw_cols
    }
    if (!is.null(loss)) {
        loss_upper <- toupper(loss)
        if (ratio) {
            metric_cols <- metric_cols[sub("_Ratio$", "", metric_cols) %in% loss_upper]
        } else {
            metric_cols <- intersect(loss_upper, metric_cols)
        }
    }
    if (length(metric_cols) == 0) stop("No matching metrics found.")

    if (type == "bar") {
        st_bar <- st
        if (ratio && is.null(methods)) {
            st_bar <- st_bar[st_bar$Method != benchmark, , drop = FALSE]
        }

        df_bar <- tidyr::pivot_longer(
            st_bar[, c("Method", metric_cols), drop = FALSE],
            cols = dplyr::all_of(metric_cols),
            names_to = "Metric", values_to = "Value"
        )

        p <- ggplot(df_bar, aes(x = .data[["Method"]], y = .data[["Value"]],
                                fill = .data[["Metric"]])) +
            geom_col(position = "dodge", width = 0.7) +
            theme_lasforecast() +
            labs(x = NULL, fill = NULL)

        if (ratio) {
            p <- p +
                geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
                labs(y = "Ratio to Benchmark")
        } else {
            p <- p + labs(y = "Value")
        }
        return(p)

    } else if (type == "heatmap") {
        df_heat <- tidyr::pivot_longer(
            st[, c("Method", metric_cols), drop = FALSE],
            cols = dplyr::all_of(metric_cols),
            names_to = "Metric", values_to = "Value"
        )
        df_heat$cell_label <- format(round(df_heat$Value, 4), nsmall = 4)

        p <- ggplot(df_heat, aes(x = .data[["Metric"]], y = .data[["Method"]],
                                 fill = .data[["Value"]])) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = .data[["cell_label"]]), size = 3) +
            scale_fill_gradient(low = "#00BFC4", high = "#F8766D", guide = "none") +
            theme_lasforecast() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
            labs(x = NULL, y = NULL, fill = NULL)
        return(p)
    }
}


#' Plot variable importance from rolling window results
#'
#' Visualizes which variables are selected across rolling windows and how
#' their coefficients evolve over time. Works with sparse methods (Lasso,
#' ALasso, TALasso, etc.) from a \code{lasforecast_roll} object.
#'
#' @param object A \code{lasforecast_roll} object.
#' @param methods Character vector of methods to visualize (default: all sparse methods).
#' @param top_n Number of top variables to show (by selection frequency). Default 20.
#'   Ignored when \code{variables} is specified.
#' @param type Type of plot: \code{"frequency"}, \code{"heatmap"}, or \code{"evolution"}.
#' @param variables Character vector of variable names to display. When
#'   specified, overrides \code{top_n}.
#' @param standardize Logical. If \code{TRUE} (default for \code{"heatmap"} and
#'   \code{"evolution"}), multiply each coefficient \eqn{\hat\beta_j} by the
#'   sample standard deviation \eqn{\hat\sigma_j} of the corresponding
#'   predictor. The resulting value represents the predicted change in \eqn{y}
#'   for a one-standard-deviation change in \eqn{x_j}, putting all coefficients
#'   on a comparable scale regardless of the variable's units. Without
#'   standardization, a few high-variance predictors can dominate the color
#'   range. Requires \code{object$x} to be available.
#' @param date_range Length-2 vector of dates (or strings coercible to Date)
#'   to restrict the plotted time span. Used by \code{"heatmap"} and
#'   \code{"evolution"} types. Ignored when dates are unavailable.
#'
#' @return A ggplot object.
#'
#' @import ggplot2
#' @export
plot_variable_importance <- function(object, methods = NULL, top_n = 20,
                                     type = c("frequency", "heatmap", "evolution"),
                                     standardize = NULL, date_range = NULL,
                                     variables = NULL) {
    type <- match.arg(type)

    if (!inherits(object, "lasforecast_roll")) {
        stop("object must be a lasforecast_roll object.")
    }

    sparse_methods <- c(
        "PLasso", "SLasso", "ALasso", "TALasso",
        "post_PLasso", "post_ALasso", "post_TALasso", "BSS"
    )
    avail <- object$methods_use
    if (is.null(methods)) {
        methods <- intersect(avail, sparse_methods)
    } else {
        methods <- intersect(methods, avail)
    }
    if (length(methods) == 0) stop("No sparse methods found in object.")

    dates <- object$dates
    has_dates <- !is.null(dates)

    # Default: standardize for heatmap/evolution, not for frequency
    if (is.null(standardize)) {
        standardize <- type %in% c("heatmap", "evolution")
    }
    # Compute variable s.d. for standardization
    x_sd <- NULL
    if (standardize) {
        if (is.null(object$x)) {
            warning("x not stored in roll object; skipping standardization.")
            standardize <- FALSE
        } else {
            x_sd <- apply(object$x, 2, sd, na.rm = TRUE)
        }
    }

    # Helper: resolve variable set from top_n or user-specified variables
    resolve_top_vars <- function(freq, top_n, variables) {
        if (!is.null(variables)) return(intersect(variables, names(freq)))
        names(sort(freq, decreasing = TRUE))[seq_len(min(top_n, length(freq)))]
    }

    if (type == "frequency") {
        freq_list <- list()
        for (m in methods) {
            beta <- object[[m]]$beta_hat
            freq <- colMeans(beta != 0)
            freq_list[[m]] <- data.frame(
                Variable = names(freq),
                Frequency = as.numeric(freq),
                Method = m,
                stringsAsFactors = FALSE
            )
        }
        df_freq <- do.call(rbind, freq_list)

        if (!is.null(variables)) {
            top_vars_union <- intersect(variables, unique(df_freq$Variable))
        } else {
            top_vars_union <- character(0)
            for (m in methods) {
                m_freq <- freq_list[[m]]
                m_sorted <- m_freq$Variable[order(m_freq$Frequency, decreasing = TRUE)]
                top_vars_union <- union(top_vars_union, m_sorted[seq_len(min(top_n, length(m_sorted)))])
            }
        }
        df_freq <- df_freq[df_freq$Variable %in% top_vars_union, ]

        # Sort by first method's frequency (descending)
        first_freq <- freq_list[[methods[1]]]
        first_freq <- first_freq[first_freq$Variable %in% top_vars_union, ]
        var_order <- first_freq$Variable[order(first_freq$Frequency, decreasing = TRUE)]
        df_freq$Variable <- factor(df_freq$Variable, levels = rev(var_order))

        p <- ggplot(df_freq, aes(x = .data[["Variable"]], y = .data[["Frequency"]], fill = .data[["Method"]])) +
            geom_col(position = "dodge", width = 0.7) +
            coord_flip() +
            theme_lasforecast() +
            labs(
                title = "Variable Selection Frequency",
                subtitle = paste0("Union of top ", top_n, " per method (", length(top_vars_union), " variables)"),
                x = NULL, y = "Selection Frequency"
            ) +
            scale_y_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1))
        return(p)
    } else if (type == "heatmap") {
        m <- methods[1]
        beta <- object[[m]]$beta_hat

        if (standardize) {
            beta <- t(t(beta) * x_sd[colnames(beta)])
        }

        # Filter by date range before selecting top vars
        time_idx <- seq_len(nrow(beta))
        if (!is.null(date_range) && has_dates) {
            date_range <- as.Date(date_range)
            time_idx <- which(dates >= date_range[1] & dates <= date_range[2])
            beta <- beta[time_idx, , drop = FALSE]
        }

        freq <- colMeans(beta != 0)
        top_vars <- resolve_top_vars(freq, top_n, variables)
        beta_sub <- beta[, top_vars, drop = FALSE]

        if (has_dates) {
            df_heat <- data.frame(Time = rep(dates[time_idx], ncol(beta_sub)))
        } else {
            df_heat <- data.frame(Time = rep(time_idx, ncol(beta_sub)))
        }
        df_heat$Variable <- rep(top_vars, each = nrow(beta_sub))
        df_heat$Value <- as.numeric(beta_sub)
        df_heat$Variable <- factor(df_heat$Variable, levels = rev(top_vars))

        fill_label <- if (standardize) "Std. Coefficient" else "Coefficient"
        p <- ggplot(df_heat, aes(x = .data[["Time"]], y = .data[["Variable"]], fill = .data[["Value"]])) +
            geom_tile() +
            scale_fill_gradient2(low = "#2196F3", mid = "white", high = "#F44336",
                                midpoint = 0, n.breaks = 5,
                                guide = guide_colorbar(barwidth = 10, barheight = 0.5)) +
            theme_lasforecast() +
            theme(axis.text.y = element_text(size = 8)) +
            labs(
                title = paste("Coefficient Heatmap -", m),
                x = if (has_dates) "Date" else "Rolling Window Index",
                y = NULL,
                fill = fill_label
            )
        return(p)
    } else if (type == "evolution") {
        time_idx <- seq_len(nrow(object[[methods[1]]]$beta_hat))
        if (!is.null(date_range) && has_dates) {
            date_range <- as.Date(date_range)
            time_idx <- which(dates >= date_range[1] & dates <= date_range[2])
        }

        beta1 <- object[[methods[1]]]$beta_hat[time_idx, , drop = FALSE]
        if (standardize) beta1 <- t(t(beta1) * x_sd[colnames(beta1)])
        freq <- colMeans(beta1 != 0)
        top_vars <- resolve_top_vars(freq, top_n, variables)

        time_vals <- if (has_dates) dates[time_idx] else time_idx
        df_list <- list()
        for (m in methods) {
            beta <- object[[m]]$beta_hat[time_idx, top_vars, drop = FALSE]
            if (standardize) beta <- t(t(beta) * x_sd[top_vars])
            df_m <- tidyr::pivot_longer(
                cbind(data.frame(Time = time_vals), as.data.frame(beta)),
                cols = -"Time", names_to = "Variable", values_to = "Coefficient"
            )
            df_m$Method <- m
            df_list[[m]] <- df_m
        }
        df_all <- do.call(rbind, df_list)
        df_all$Variable <- factor(df_all$Variable, levels = top_vars)

        y_label <- if (standardize) "Std. Coefficient" else "Coefficient"
        p <- ggplot(df_all, aes(x = .data[["Time"]], y = .data[["Coefficient"]],
                                color = .data[["Method"]])) +
            geom_line(alpha = 0.7) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
            facet_wrap(~ Variable, scales = "free_y") +
            theme_lasforecast() +
            labs(
                title = "Coefficient Path",
                x = if (has_dates) "Date" else "Rolling Window Index",
                y = y_label
            )
        return(p)
    }
}
