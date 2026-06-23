#' Print method for lasforecast_model objects
#'
#' @param x A \code{lasforecast_model} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.lasforecast_model <- function(x, ...) {
    cat("LasForecast Model\n")
    cat("  Method:", x$method, "\n")

    if (identical(x$method, "arma")) {
        cat("  ARMA order: (", paste(x$arma_order, collapse = ", "), ")\n", sep = "")
        cat("  Parameters:", x$df, "\n")
        cat("  Forecast:", format(x$forecast, digits = 4), "\n")
        return(invisible(x))
    }

    coefs <- x$coefficients
    n_coef <- length(coefs)
    n_nonzero <- sum(coefs != 0)

    cat("  Coefficients:", n_coef, "(", n_nonzero, "nonzero )\n")

    if (!is.null(x$lambda)) {
        cat("  Lambda:", format(x$lambda, digits = 4), "\n")
    }
    if (!is.null(x$tuning_param)) {
        cat("  Tuning param:", x$tuning_param, "\n")
    }

    invisible(x)
}


#' Extract coefficients from a lasforecast_model
#'
#' @param object A \code{lasforecast_model} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Numeric vector of coefficients.
#' @export
coef.lasforecast_model <- function(object, ...) {
    object$coefficients
}


#' Print method for lasforecast_roll objects
#'
#' @param x A \code{lasforecast_roll} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.lasforecast_roll <- function(x, ...) {
    methods <- x$methods_use
    num_forecast <- length(x[[methods[1]]]$y_hat)

    cat("LasForecast Rolling Window Forecast\n")
    cat("  Methods:", paste(methods, collapse = ", "), "\n")
    cat("  Window size:", x$roll_window, "\n")
    if (!is.null(x$window_type)) {
        cat("  Window type:", x$window_type, "\n")
    }
    cat("  Horizon:", x$h, "\n")
    cat("  Forecast periods:", num_forecast, "\n")

    if (!is.null(x$dates)) {
        cat("  Date range:", format(min(x$dates)), "to", format(max(x$dates)), "\n")
    }

    cat("\n  MSE:\n")
    mse_df <- data.frame(
        Method = names(x$mse),
        MSE = as.numeric(x$mse),
        RMSE = sqrt(as.numeric(x$mse))
    )
    print(mse_df, row.names = FALSE, digits = 4)

    invisible(x)
}


#' Summary method for lasforecast_roll objects
#'
#' @param object A \code{lasforecast_roll} object.
#' @param benchmark Benchmark method for ratio computation (default \code{"RW"}).
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame with performance summary.
#' @export
summary.lasforecast_roll <- function(object, benchmark = "RW", ...) {
    compare_forecasts(object, benchmark = benchmark)
}


#' Print method for lasforecast_backtest objects
#'
#' @param x A \code{lasforecast_backtest} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.lasforecast_backtest <- function(x, ...) {
    cat("LasForecast Backtest\n")
    cat("  Methods:", paste(x$methods, collapse = ", "), "\n")
    cat("  Benchmark:", x$benchmark, "\n")
    cat("  Horizon:", x$h, "\n")
    cat("  Window:", x$roll_window, "(", x$window_type, ")\n")
    cat("  Loss functions:", paste(x$loss, collapse = ", "), "\n")

    if (!is.null(x$dates)) {
        cat("  Date range:", format(min(x$dates)), "to", format(max(x$dates)), "\n")
    }

    cat("\nPerformance Summary:\n")
    print(x$summary_table, row.names = FALSE, digits = 4)

    invisible(x)
}


#' Summary method for lasforecast_backtest objects
#'
#' @param object A \code{lasforecast_backtest} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame with performance summary including benchmark ratios.
#' @export
summary.lasforecast_backtest <- function(object, ...) {
    st <- object$summary_table
    bench <- object$benchmark
    loss_cols <- setdiff(names(st), c("Method", grep("_Ratio$", names(st), value = TRUE)))

    for (col in loss_cols) {
        bench_val <- st[[col]][st$Method == bench]
        if (length(bench_val) == 1 && bench_val != 0) {
            st[[paste0(col, "_Ratio")]] <- st[[col]] / bench_val
        }
    }

    st
}


# --- lasforecast_debias_ivx methods ---

#' Print method for lasforecast_debias_ivx objects
#'
#' @param x A \code{lasforecast_debias_ivx} object.
#' @param ... Additional arguments (ignored).
#'
#' @export
print.lasforecast_debias_ivx <- function(x, ...) {
    b_las <- as.numeric(x$b_hat_las)

    cat("IVX-Desparsified Lasso Inference\n")
    cat("  Focal variables:", paste(x$term_names[x$d_ind], collapse = ", "), "\n")
    cat("  Lasso nonzero:", sum(b_las != 0), "of", length(b_las), "\n")
    cat("  Lambda (first step):", format(x$lambda_hat[1], digits = 4), "\n")

    invisible(x)
}


#' Extract focal IVX coefficients from a lasforecast_debias_ivx object
#'
#' @param object A \code{lasforecast_debias_ivx} object.
#' @param ... Additional arguments (ignored).
#'
#' @return Named numeric vector of debiased IVX estimates for the focal variables.
#' @export
coef.lasforecast_debias_ivx <- function(object, ...) {
    est <- object$theta_hat_ivx
    names(est) <- object$term_names[object$d_ind]
    est
}


#' Summary method for lasforecast_debias_ivx objects
#'
#' Produces an \code{lm}-style coefficient table for the focal variables
#' (debiased IVX estimates with standard errors, t values, and p-values),
#' plus the nonzero first-step Lasso coefficients for control variables.
#'
#' @param object A \code{lasforecast_debias_ivx} object.
#' @param ... Additional arguments (ignored).
#'
#' @return An object of class \code{"summary.lasforecast_debias_ivx"}.
#'
#' @importFrom stats pnorm printCoefmat
#' @export
summary.lasforecast_debias_ivx <- function(object, ...) {
    focal_names <- object$term_names[object$d_ind]

    est <- object$theta_hat_ivx
    se <- object$sigma_hat_ivx
    tval <- est / se
    pval <- 2 * pnorm(abs(tval), lower.tail = FALSE)

    focal_mat <- cbind(
        Estimate = est, `Std. Error` = se,
        `t value` = tval, `Pr(>|t|)` = pval
    )
    rownames(focal_mat) <- focal_names

    # Nonzero Lasso controls (excluding focal variables)
    b_las <- as.numeric(object$b_hat_las)
    names(b_las) <- object$term_names
    nz_control_ind <- setdiff(which(b_las != 0), object$d_ind)
    lasso_controls <- b_las[nz_control_ind]

    structure(
        list(
            focal = focal_mat,
            lasso_controls = lasso_controls,
            wald_stat = object$wald_stat,
            p_value_wald = object$p_value_wald,
            lambda_hat = object$lambda_hat
        ),
        class = "summary.lasforecast_debias_ivx"
    )
}


#' Print method for summary.lasforecast_debias_ivx objects
#'
#' @param x A \code{summary.lasforecast_debias_ivx} object.
#' @param ... Additional arguments passed to \code{\link[stats]{printCoefmat}}.
#'
#' @export
print.summary.lasforecast_debias_ivx <- function(x, ...) {
    cat("IVX-Desparsified Lasso Inference\n---\n")

    cat("Focal coefficients:\n")
    printCoefmat(x$focal, signif.stars = TRUE, has.Pvalue = TRUE, ...)

    if (!is.na(x$wald_stat)) {
        cat("\nWald test: chi2(", nrow(x$focal), ") = ",
            format(x$wald_stat, digits = 4), ", p = ",
            format(x$p_value_wald, digits = 4), "\n", sep = "")
    }

    if (length(x$lasso_controls) > 0) {
        cat("\nLasso-selected controls (first-step estimates):\n")
        print(round(x$lasso_controls, 6))
    }

    invisible(x)
}
