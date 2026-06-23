#' Rolling window forecast (Refactored)
#'
#' The function reads data and make forecasts based on linear predictive regression
#' with diverse methods. It incorporates both short-horizon and long-horizon forecasting.
#'
#' @param x Full sample predictor. Can be a matrix, data.frame, tibble, or tsibble. Should not include the intercept column.
#' @param y Full sample forecast target, or column name string when \code{x} is a data frame.
#'   Note: in the code, y(t+1) is predicted by x(t). If you want to do long horizon forecast (h > 1), you must manually shift y accordingly before input,
#'   i.e. before input, change y(t) to y(t+h-1), such that in the code y(t+h) is predicted by x(t).
#' @param roll_window Length of the rolling window
#' @param h Forecast horizon. The rolling regression is always one-step (it pairs
#'   \code{y(t+1)} with \code{x(t)}); \code{h} only sets the gap between the last
#'   training target \code{y(i-h)} and the forecast input \code{x(i-1)}, which
#'   prevents look-ahead when \code{y} is (or has been shifted to be) an
#'   \code{h}-horizon target. See the \code{y} note on shifting.
#' @param methods_use method in use
#' @param train_method parameter tuning method for Lasso type methods.
#' @param verbose boolean to control whether print information on screen
#' @param ar_order 0 or 1 to control whether include ar1 lag or not
#' @param date_col Optional date column name (for data.frame/tibble input).
#' @param window_type Character string: \code{"rolling"} (default) for fixed-size
#'   rolling window, or \code{"expanding"} for expanding window.
#' @param type Type of forecast target y. This is only relevant for RW method.
#'  "diff" for returns/differences and RW predicts zeros; "level" for levels and RW predicts the previous value. Default is "level".
#' @param alasso_model Initial estimation model for adaptive Lasso: \code{"lasso"}, \code{"ridge"}, or \code{"ols"}. Default \code{NULL} auto-selects based on dimensionality.
#' @param alasso_lambda_init Penalty for the Lasso/ridge initial estimator
#'   generating the adaptive Lasso weights (ALasso/TALasso and their post
#'   variants). Default \code{NULL}: tuned internally by 10-fold block CV
#'   within each window. Ignored when the initial estimator is OLS.
#' @param arma_order ARMA order \code{c(p, d, q)} for the \code{"ARMA"} benchmark.
#'   Default \code{NULL} selects \eqn{p} and \eqn{q} automatically by information
#'   criterion. ARMA is univariate (it ignores the predictors and \code{ar_order});
#'   other ARMA specifications (\code{include.mean}, \code{d}, \code{max_p},
#'   \code{max_q}, \code{ic}, \code{arma_method}) may be passed via \code{...}.
#'   See \code{\link{arma}}.
#' @param bss_k Fixed subset size for the \code{"BSS"} method. Default
#'   \code{NULL}: the subset size is selected within each window by
#'   information criterion (\code{train_method} when it is one of
#'   \code{"aic"}, \code{"bic"}, \code{"aicc"}, \code{"hqc"}; otherwise BIC).
#'   \code{"BSS"} requires the Gurobi solver; see \code{\link{bss}}, whose
#'   other options (\code{ic_type}, \code{k_max}, \code{time_limit}) may be
#'   passed via \code{...}.
#' @param ... Additional arguments passed to \code{\link{lasso}} and its
#'   tuning function \code{train_lasso}. Useful options include
#'   \code{intercept}, \code{gamma}, \code{nlambda}, \code{k},
#'   \code{alpha}, \code{lambda_seq}, \code{lambda_min_ratio}.
#'   Note \code{k} is the fold count for cross-validation tuning of the
#'   Lasso methods, not the BSS subset size (\code{bss_k}).
#'
#' @import glmnet
#' @export
#'
roll_predict <- function(x, y = NULL, roll_window, h = 1,
                         methods_use = c(
                             "RW", "RWwD", "OLS", "PLasso", "SLasso",
                             "ALasso", "TALasso"
                         ),
                         train_method = "cv", verbose = TRUE, ar_order = 0,
                         date_col = NULL, window_type = c("rolling", "expanding"), type = "level",
                         alasso_model = NULL, alasso_lambda_init = NULL, arma_order = NULL,
                         bss_k = NULL, ...) {
    window_type <- match.arg(window_type)

    inp <- prepare_input(x, y, date_col)
    x <- inp$x_mat
    y <- inp$y_vec
    dates_all <- inp$dates

    methods_all <- c("RW", "RWwD", "OLS", "ARMA", "PLasso", "SLasso", "ALasso", "TALasso",
                      "post_PLasso", "post_SLasso", "post_ALasso", "post_TALasso", "BSS")
    invalid_methods <- setdiff(methods_use, methods_all)
    if (length(invalid_methods) > 0) {
        stop(
            "Invalid methods in methods_use: ",
            paste(invalid_methods, collapse = ", "),
            ". Method names are case-sensitive."
        )
    }
    if ("BSS" %in% methods_use && !requireNamespace("gurobi", quietly = TRUE)) {
        stop("Method 'BSS' requires the 'gurobi' package (shipped with the ",
             "Gurobi Optimizer, https://www.gurobi.com).")
    }

    # Arguments pinned inside the per-method dispatch must not also arrive via
    # `...`: the duplicate-argument error would be demoted to a per-window
    # warning by the tryCatch below, and the method would silently return NA
    # forecasts instead of failing.
    dots_names <- names(list(...))
    if ("model" %in% dots_names) {
        stop("Argument 'model' is reserved; use 'alasso_model' to set the adaptive Lasso initial estimator.")
    }
    if ("lambda_init" %in% dots_names) {
        stop("Argument 'lambda_init' is reserved; use 'alasso_lambda_init' to set the initial estimator's penalty.")
    }
    if ("predictive" %in% dots_names) {
        stop("Argument 'predictive' is set internally; the rolling loop already aligns y and x.")
    }
    scale_pinned <- c("PLasso", "SLasso", "post_PLasso", "post_SLasso")
    if ("scale_x" %in% dots_names && any(methods_use %in% scale_pinned)) {
        stop("Argument 'scale_x' cannot be combined with method(s) ",
             paste(intersect(methods_use, scale_pinned), collapse = ", "),
             ", whose scale convention is fixed by definition. Drop 'scale_x' or those methods.")
    }

    n <- nrow(x)
    p <- ncol(x)
    num_methods <- length(methods_use)

    if ("OLS" %in% methods_use && p > roll_window) {
        stop("Error: OLS is not feasible when p > roll_window. Please remove 'OLS' from methods_use or reduce the number of predictors.")
    }

    # Forecast count
    if (ar_order == 0) {
        num_forecast <- n - roll_window
    } else {
        num_forecast <- n - roll_window - 2 * h + 1
        if (h == 1) {
            warning("Warning: AR term with h=1 might cause multicollinearity with dividend yields.")
        }
    }

    # Initialize Result Containers
    # We use a list of lists structure
    save_result <- vector("list", num_methods)
    names(save_result) <- methods_use

    # Define column names for beta matrix
    if (ar_order == 0) {
        beta_cols <- colnames(x)
    } else {
        beta_cols <- c(colnames(x), "AR(1)")
    }
    p_eff <- p + ar_order

    for (m in methods_use) {
        save_result[[m]] <- list(
            # NA (not 0) so a window where the method errors out is dropped from
            # the loss metrics rather than scored as if it forecast exactly 0
            # (a real 0 forecast, e.g. RW on differences, overwrites this).
            y_hat = rep(NA_real_, num_forecast), # forecasts
            beta_hat = matrix(0, num_forecast, p_eff, dimnames = list(NULL, beta_cols)),
            tuning_param = rep(NA, num_forecast), # Store lambda/k etc if applicable
            df = rep(NA, num_forecast) # degrees of freedom
        )
    }

    # Prediction Loop
    t0 <- roll_window + 1

    for (i in t0:n) {
        # Progress
        # With an AR term the first feasible origin is i = roll_window + 2h
        # (the y-lag must reach back to index 1), so the slot index shifts by
        # 2h - 1, not 2h: the first processed window must land in slot 1.
        tt <- i - roll_window - ar_order * (2 * h - 1)
        if (verbose && (tt %% 10 == 0 || tt == 1)) {
            cat("Rolling Window =", roll_window, ", h =", h, ", Prediction:", tt, "/", num_forecast, "\n")
        }

        t_start <- Sys.time()

        # --- Data Preparation ---
        if (ar_order == 0) {
            # The regression always pairs y[t+1] ~ x[t] (one-step in the input
            # series). The h-period gap -- training targets end at y[i-h] while
            # the forecast uses x[i-1] -- is what makes this a leakage-free
            # h-step forecast: it presumes y already realises h-1 periods after
            # its index (a genuine long-horizon target, or y shifted by
            # backtest(shift_y = TRUE)). The horizon is not applied twice.
            if (i < t0 + h || window_type == "expanding") {
                x_est <- as.matrix(x[1:(i - h - 1), ])
                y_est <- as.matrix(y[2:(i - h)])
                x_for <- x[i - 1, ]
            } else {
                x_est <- as.matrix(x[(i - roll_window - h):(i - h - 1), ])
                y_est <- as.matrix(y[(i - roll_window - h + 1):(i - h)])
                x_for <- x[i - 1, ]
            }
        } else {
            # AR(1) logic
            if (i < t0 + 2 * h - 1) {
                next
            } else if (window_type == "expanding") {
                # Pair y_t with x_{t-1} and the lag y_{t-h}; training rows
                # start at t = h + 1 so the h-lag exists for every row
                # (for h = 1 this reduces to the previous indexing).
                x_est <- cbind(
                    as.matrix(x[h:(i - h - 1), ]),
                    as.matrix(y[1:(i - 2 * h)])
                )
                y_est <- as.matrix(y[(h + 1):(i - h)])
                x_for <- c(x[i - 1, ], y[i - h])
            } else {
                # Add lagged y to x
                # x_est includes lags of x AND lags of y
                x_est <- cbind(
                    as.matrix(x[(i - roll_window - h):(i - h - 1), ]),
                    as.matrix(y[(i - roll_window - 2 * h + 1):(i - 2 * h)])
                )
                y_est <- as.matrix(y[(i - roll_window - h + 1):(i - h)])
                x_for <- c(x[i - 1, ], y[i - h])
            }
        }

        # --- Estimation & Prediction ---

        for (method in methods_use) {
            # Dispatch to specific estimate function
            # We map method string to function call

            model <- tryCatch(
                {
                    switch(method,
                        "RW" = rw(x_est, y_est, type = type),
                        "RWwD" = rwwd(x_est, y_est),
                        "OLS" = ols(x_est, y_est),
                        "ARMA" = {
                            # Univariate ARMA on the response window; other specs
                            # (include.mean, d, max_p/q, ic, arma_method) flow via ...
                            arma_dots <- list(...)
                            arma_dots$order <- NULL # reserved for arma_order
                            do.call(arma, c(list(y = y_est, order = arma_order, h = h), arma_dots))
                        },
                        "PLasso" = lasso(x_est, y_est, method = "lasso", predictive = FALSE, scale_x = FALSE, train_method = train_method, ...),
                        "SLasso" = lasso(x_est, y_est, method = "lasso", predictive = FALSE, scale_x = TRUE, train_method = train_method, ...),
                        "ALasso" = lasso(x_est, y_est, method = "alasso", predictive = FALSE, train_method = train_method, model = alasso_model, lambda_init = alasso_lambda_init, ...),
                        "TALasso" = lasso(x_est, y_est, method = "talasso", predictive = FALSE, train_method = train_method, model = alasso_model, lambda_init = alasso_lambda_init, ...),
                        "post_PLasso" = lasso(x_est, y_est, method = "post_lasso", predictive = FALSE, scale_x = FALSE, train_method = train_method, ...),
                        "post_SLasso" = lasso(x_est, y_est, method = "post_lasso", predictive = FALSE, scale_x = TRUE, train_method = train_method, ...),
                        "post_ALasso" = lasso(x_est, y_est, method = "post_alasso", predictive = FALSE, train_method = train_method, model = alasso_model, lambda_init = alasso_lambda_init, ...),
                        "post_TALasso" = lasso(x_est, y_est, method = "post_talasso", predictive = FALSE, train_method = train_method, model = alasso_model, lambda_init = alasso_lambda_init, ...),
                        "BSS" = {
                            # Subset size: bss_k if fixed, else per-window IC
                            # selection driven by train_method when it is an
                            # IC. 'k' in `...` is the CV fold count for the
                            # Lasso methods, not the subset size: drop it.
                            bss_dots <- list(...)
                            bss_dots$k <- NULL
                            if (is.null(bss_dots$ic_type) &&
                                    train_method %in% c("aic", "bic", "aicc", "hqc")) {
                                bss_dots$ic_type <- train_method
                            }
                            do.call(bss, c(
                                list(x = x_est, y = y_est, k = bss_k, predictive = FALSE),
                                bss_dots
                            ))
                        },
                        stop(paste("Unknown method:", method))
                    )
                },
                error = function(e) {
                    warning(paste("Error in method", method, ":", e$message))
                    return(NULL)
                }
            )

            if (!is.null(model)) {
                if (identical(model$method, "arma")) {
                    # Univariate benchmark: no predictor coefficients to store.
                    save_result[[method]]$beta_hat[tt, ] <- NA_real_
                    save_result[[method]]$df[tt] <- model$df
                } else {
                    # Store fit details
                    # Helper to safely extract elements
                    safe_get <- function(obj, name, default = NA) if (!is.null(obj[[name]])) obj[[name]] else default

                    # Coefficients
                    coefs <- model$coefficients
                    # Beta hat usually excludes intercept in save_result structure?
                    # Original code: save_result$OLS$beta_hat[tt, ] <- coef_ols[-1]
                    # Most methods return intercept as first element.
                    # But bss returns intercept inside?
                    # Let's assume standard is Intercept First.
                    if (length(coefs) == p_eff + 1) {
                        beta_hat <- coefs[-1]
                    } else if (length(coefs) == p_eff) {
                        beta_hat <- coefs
                    } else {
                        warning("Unexpected coefficient length for method ", method)
                        beta_hat <- rep(NA, p_eff)
                    }

                    save_result[[method]]$beta_hat[tt, ] <- beta_hat

                    # Tuning Param
                    if (!is.null(model$lambda)) save_result[[method]]$tuning_param[tt] <- model$lambda
                    if (!is.null(model$tuning_param)) save_result[[method]]$tuning_param[tt] <- model$tuning_param # for bss

                    # DF
                    # df is number of non-zero slopes
                    save_result[[method]]$df[tt] <- sum(beta_hat != 0)
                }

                # Prediction
                y_pred <- predict(model, newx = x_for)
                save_result[[method]]$y_hat[tt] <- y_pred

            }
        }

        t_end <- Sys.time()
        # if(verb) print(t_end - t_start)
    }

    # MSE Calculation
    y_0 <- y[-(1:roll_window)]
    # Adjust y_0 length if needed?
    # Original code: y_0 <- y[-(1:roll_window)]. Length n - roll_window.
    # If ar_order=1, we skip more?
    # num_forecast is calculated.
    # If ar_order=1, tt index matches.

    # For safety, let's use the valid indices.
    # But y_hat is filled.
    # We should align y_0.

    # Replicating original mse logic:
    mse <- rep(0, num_methods)
    names(mse) <- methods_use

    # We need the actual values corresponding to the forecasts.
    # Logic: y_hat[tt] corresponds to prediction for y[i]?
    # Original: y_est = y[...]; x_for = x[i-1]; predict for y[i].
    # So y_hat[tt] is forecast for y[i].

    # y_target vector construction
    y_target <- numeric(num_forecast)
    for (i in t0:n) {
        tt <- i - roll_window - ar_order * (2 * h - 1)
        if (tt > 0 && tt <= num_forecast) {
            y_target[tt] <- y[i]
        }
    }

    for (j in methods_use) {
        mse[j] <- mean((save_result[[j]]$y_hat - y_target)^2, na.rm = TRUE)
    }

    save_result$mse <- mse
    save_result$y <- y
    save_result$x <- x
    save_result$methods_use <- methods_use
    save_result$roll_window <- roll_window
    save_result$h <- h
    save_result$window_type <- window_type

    # Store dates aligned with forecasts
    if (!is.null(dates_all) && length(dates_all) == n) {
        save_result$dates <- tail(dates_all, num_forecast)
    }

    class(save_result) <- "lasforecast_roll"
    return(save_result)
}
