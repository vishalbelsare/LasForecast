#' Backtest Forecasting Methods
#'
#' Runs rolling window forecasts for multiple methods and compares their
#' out-of-sample prediction performance. This is the main entry point for
#' benchmarking in LasForecast.
#'
#' @param x Predictor data. Can be a matrix, data.frame, tibble, or tsibble.
#' @param y Response vector, or column name string when \code{x} is a data frame.
#' @param methods Character vector of methods to compare. Default includes
#'   common methods: \code{"RW"}, \code{"OLS"}, \code{"PLasso"}, \code{"ALasso"}.
#'   See \code{\link{roll_predict}} for all supported methods.
#' @param roll_window Rolling window size. If NULL, defaults to floor(n * 0.5).
#' @param h Forecast horizon (default 1).
#' @param window_type \code{"rolling"} (default) or \code{"expanding"}.
#' @param loss Character vector of loss functions to compute. Options:
#'   \code{"mse"}, \code{"rmse"}, \code{"mae"}, \code{"mape"}, \code{"r2oos"}.
#'   \code{"r2oos"} is the Campbell-Thompson out-of-sample \eqn{R^2}, computed
#'   relative to \code{benchmark} (\eqn{1 - \sum e_{method}^2 / \sum e_{bench}^2}).
#' @param benchmark Benchmark method name (default \code{"RWwD"}). Used as the
#'   reference forecast for both \code{RMSE_Ratio} and the \code{"r2oos"} metric.
#' @param train_method Training method for Lasso-type parameter tuning (default \code{"cv"}).
#' @param ar_order 0 or 1, whether to include AR(1) lag.
#' @param date_col Optional date column name.
#' @param verbose Logical, whether to print progress (default TRUE).
#' @param type Type of forecast target y. This is only relevant for RW method.
#'  "diff" for returns/differences and RW predicts zeros; "level" for levels and RW predicts the previous value. Default is "level".
#' @param alasso_model Initial estimation model for adaptive Lasso: \code{"lasso"}, \code{"ridge"}, or \code{"ols"}. Default \code{NULL} auto-selects based on dimensionality.
#' @param alasso_lambda_init Penalty for the Lasso/ridge initial estimator
#'   generating the adaptive Lasso weights (ALasso/TALasso and their post
#'   variants). Default \code{NULL}: tuned internally by 10-fold block CV
#'   within each window. Ignored when the initial estimator is OLS.
#' @param arma_order ARMA order \code{c(p, d, q)} for the \code{"ARMA"} benchmark.
#'   Default \code{NULL} selects \eqn{p} and \eqn{q} automatically by information
#'   criterion. Other ARMA specifications (\code{include.mean}, \code{d},
#'   \code{max_p}, \code{max_q}, \code{ic}, \code{arma_method}) may be passed via
#'   \code{...}. See \code{\link{arma}}.
#' @param bss_k Fixed subset size for the \code{"BSS"} method. Default
#'   \code{NULL}: selected within each window by information criterion
#'   (\code{train_method} when it is an IC name, otherwise BIC).
#'   \code{"BSS"} requires the Gurobi solver; see \code{\link{bss}}.
#' @param shift_y Logical (default FALSE). Controls how the horizon \code{h} is
#'   applied when \code{h > 1}:
#'   \describe{
#'     \item{\code{TRUE}}{\code{y} is an ordinary one-period series and is
#'       shifted forward by \code{h - 1} so the rolling regression targets
#'       \code{y(t + h)} with \code{x(t)} (direct multi-step). Use this to obtain
#'       an \code{h}-step-ahead forecast of a one-period series.}
#'     \item{\code{FALSE}}{\code{y} is assumed to be \emph{already} an
#'       \code{h}-horizon target that realises at time \code{t + h - 1}
#'       (e.g., \code{y_t = log P_{t+h-1} - log P_{t-1}}); the rolling loop's
#'       internal \code{h}-period gap then guards against look-ahead. Leaving
#'       the default with a one-period \code{y} and \code{h > 1} does \emph{not}
#'       give an \code{h}-step forecast.}
#'   }
#'   The horizon is never applied twice: the shift fixes the regression target
#'   and the internal gap prevents leakage, two distinct roles.
#' @param ... Additional arguments passed to \code{\link{lasso}} via
#'   \code{\link{roll_predict}}. Useful options include
#'   \code{intercept}, \code{gamma}, \code{nlambda}, \code{k},
#'   \code{alpha}, \code{lambda_seq}, \code{lambda_min_ratio}.
#'
#' @return An object of class \code{"lasforecast_backtest"} containing:
#' \describe{
#'   \item{results}{A \code{lasforecast_roll} object with rolling forecast results.}
#'   \item{summary_table}{Data frame with performance metrics for all methods.}
#'   \item{dates}{Date vector aligned with forecasts (if available).}
#'   \item{benchmark}{Name of the benchmark method.}
#'   \item{methods}{Methods compared.}
#'   \item{loss}{Loss functions used.}
#' }
#'
#' @examples
#' \dontrun{
#' data("fredmd")
#' df <- na.omit(fredmd)
#' bt <- backtest(df,
#'     y = "infl_cpi",
#'     methods = c("RW", "OLS", "PLasso", "ALasso"),
#'     roll_window = 120, h = 1
#' )
#' print(bt)
#' autoplot(bt)
#' }
#'
#' @import ggplot2 dplyr tidyr
#' @importFrom rlang .data
#' @importFrom dplyr all_of
#' @importFrom utils tail
#' @importFrom stats median coef predict lsfit sd rnorm cov lm pchisq toeplitz
#' @importFrom utils combn
#' @importFrom methods as
#' @export
backtest <- function(x, y = NULL,
                     methods = c("RW", "OLS", "PLasso", "ALasso"),
                     roll_window = NULL,
                     h = 1,
                     window_type = c("rolling", "expanding"),
                     loss = c("mse", "rmse", "mae"),
                     benchmark = "RWwD",
                     train_method = "cv",
                     ar_order = 0,
                     date_col = NULL,
                     verbose = TRUE,
                     type = "level",
                     alasso_model = NULL,
                     alasso_lambda_init = NULL,
                     arma_order = NULL,
                     bss_k = NULL,
                     shift_y = FALSE, ...) {
    window_type <- match.arg(window_type)

    inp <- prepare_input(x, y, date_col)
    x_mat <- inp$x_mat
    y_vec <- inp$y_vec
    dates_all <- inp$dates
    n <- nrow(x_mat)

    # Internal y-shifting for direct multi-step forecasting
    if (shift_y && h > 1) {
        n_new <- n - h + 1
        y_vec <- y_vec[h:n]
        x_mat <- x_mat[1:n_new, , drop = FALSE]
        # Align dates with forecast targets (y_orig[t+h-1]), not x rows
        if (!is.null(dates_all)) dates_all <- dates_all[h:n]
        n <- n_new
    } else if (!shift_y && h > 1) {
        warning(
            "h > 1 with shift_y = FALSE assumes 'y' is already an h-horizon ",
            "target (realising at t + h - 1). If 'y' is an ordinary one-period ",
            "series, set shift_y = TRUE to obtain a direct h-step forecast; ",
            "otherwise the rolling regression stays one-step and simply drops ",
            "the most recent h - 1 usable observations."
        )
    }

    if (is.null(roll_window)) {
        roll_window <- floor(n * 0.5)
    }

    if (!benchmark %in% methods) {
        methods <- c(benchmark, methods)
    }

    res <- roll_predict(
        x = x_mat, y = y_vec,
        roll_window = roll_window,
        h = h,
        methods_use = methods,
        train_method = train_method,
        verbose = verbose,
        ar_order = ar_order,
        window_type = window_type,
        type = type,
        alasso_model = alasso_model,
        alasso_lambda_init = alasso_lambda_init,
        arma_order = arma_order,
        bss_k = bss_k,
        ...
    )

    if (!is.null(dates_all)) {
        num_forecast <- length(res[[methods[1]]]$y_hat)
        res$dates <- tail(dates_all, num_forecast)
    }

    # Calculate metrics
    num_forecast <- length(res[[methods[1]]]$y_hat)
    y_target <- tail(y_vec, num_forecast)

    bench_pred <- res[[benchmark]]$y_hat
    summary_rows <- lapply(methods, function(m) {
        metrics <- calc_loss_all(y_target, res[[m]]$y_hat, loss,
                                 benchmark_pred = bench_pred)
        row <- data.frame(Method = m, stringsAsFactors = FALSE)
        for (lname in names(metrics)) {
            row[[lname]] <- metrics[[lname]]
        }
        row
    })
    summary_table <- do.call(rbind, summary_rows)

    if (benchmark %in% methods && "RMSE" %in% names(summary_table)) {
        bench_rmse <- summary_table$RMSE[summary_table$Method == benchmark]
        summary_table$RMSE_Ratio <- summary_table$RMSE / bench_rmse
    }

    structure(
        list(
            results = res,
            summary_table = summary_table,
            dates = res$dates,
            benchmark = benchmark,
            methods = methods,
            loss = loss,
            h = h,
            roll_window = roll_window,
            window_type = window_type
        ),
        class = "lasforecast_backtest"
    )
}


#' Calculate a single loss metric
#'
#' @param actual Numeric vector of actual values.
#' @param predicted Numeric vector of predicted values.
#' @param type Loss function type: \code{"mse"}, \code{"rmse"}, \code{"mae"},
#'   \code{"mape"}, or \code{"r2oos"}.
#' @param benchmark_pred Numeric vector of benchmark forecasts, aligned with
#'   \code{actual}. Required only for \code{"r2oos"}: it supplies the denominator
#'   of the Campbell-Thompson out-of-sample \eqn{R^2}. Returns \code{NA} for
#'   \code{"r2oos"} when not supplied.
#'
#' @return Scalar loss value.
#' @keywords internal
calc_loss <- function(actual, predicted, type = "mse", benchmark_pred = NULL) {
    e <- actual - predicted
    switch(type,
        "mse" = mean(e^2, na.rm = TRUE),
        "rmse" = sqrt(mean(e^2, na.rm = TRUE)),
        "mae" = mean(abs(e), na.rm = TRUE),
        "mape" = {
            nonzero <- actual != 0
            if (sum(nonzero) == 0) {
                return(NA_real_)
            }
            mean(abs(e[nonzero] / actual[nonzero]), na.rm = TRUE) * 100
        },
        "r2oos" = {
            # Campbell & Thompson (2008) out-of-sample R^2: one minus the ratio of
            # the method's sum of squared forecast errors to the benchmark's.
            # R2_OOS > 0 means the method forecasts better than the benchmark.
            # Numerator and denominator are summed over the same windows (those
            # where both forecasts are available) so the ratio stays a like-for-
            # like comparison when some windows are missing.
            if (is.null(benchmark_pred)) {
                return(NA_real_)
            }
            e_bench <- actual - benchmark_pred
            ok <- is.finite(e) & is.finite(e_bench)
            if (!any(ok)) {
                return(NA_real_)
            }
            ss_res <- sum(e[ok]^2)
            ss_bench <- sum(e_bench[ok]^2)
            if (ss_bench == 0) {
                return(NA_real_)
            }
            1 - ss_res / ss_bench
        },
        stop("Unknown loss type: ", type)
    )
}


#' Calculate all requested loss metrics
#'
#' @param actual Numeric vector of actual values.
#' @param predicted Numeric vector of predicted values.
#' @param loss Character vector of loss types.
#' @param benchmark_pred Numeric vector of benchmark forecasts (for
#'   \code{"r2oos"}); see \code{\link{calc_loss}}.
#'
#' @return Named list of loss values.
#' @keywords internal
calc_loss_all <- function(actual, predicted, loss = c("mse", "rmse", "mae"),
                          benchmark_pred = NULL) {
    metrics <- list()
    for (l in loss) {
        val <- calc_loss(actual, predicted, l, benchmark_pred = benchmark_pred)
        # Use uppercase names for display
        name <- toupper(l)
        if (name == "R2OOS") name <- "R2_OOS"
        metrics[[name]] <- val
    }
    metrics
}


#' Compare Forecast Performance
#'
#' Evaluates and compares the performance of different forecasting methods
#' from a \code{lasforecast_roll} object. Use \code{autoplot()} for plots.
#'
#' @param result_list Output from \code{roll_predict}.
#' @param benchmark Name of the benchmark method (default is "RW" or first method).
#'   Used as the reference forecast for both \code{RMSE_Ratio} and the
#'   \code{"r2oos"} (Campbell-Thompson) metric.
#' @param loss Character vector of loss functions. Options:
#'   \code{"mse"}, \code{"rmse"}, \code{"mae"}, \code{"mape"}, \code{"r2oos"}.
#'
#' @return A data.frame with accuracy metrics per method.
#'
#' @import ggplot2 dplyr tidyr
#' @export
compare_forecasts <- function(result_list, benchmark = "RW",
                              loss = c("rmse", "mae")) {
    methods <- result_list$methods_use
    y <- result_list$y

    num_forecast <- length(result_list[[methods[1]]]$y_hat)
    y_target <- tail(y, num_forecast)

    bench_pred <- if (benchmark %in% methods) result_list[[benchmark]]$y_hat else NULL
    summary_rows <- lapply(methods, function(m) {
        metrics <- calc_loss_all(y_target, result_list[[m]]$y_hat, loss,
                                 benchmark_pred = bench_pred)
        row <- data.frame(Method = m, stringsAsFactors = FALSE)
        for (lname in names(metrics)) {
            row[[lname]] <- metrics[[lname]]
        }
        row
    })
    summary_stats <- do.call(rbind, summary_rows)

    if (benchmark %in% methods && "RMSE" %in% names(summary_stats)) {
        bench_rmse <- summary_stats$RMSE[summary_stats$Method == benchmark]
        summary_stats$RMSE_Ratio <- summary_stats$RMSE / bench_rmse
    }

    summary_stats
}
