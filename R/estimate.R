#' Lasso-Family Estimation for Predictive Regression
#'
#' Tunes the penalty parameter and fits a Lasso-family model.
#' Supports the standard Lasso, Adaptive Lasso (ALasso), Twin Adaptive Lasso
#' (TALasso), and Post-Lasso OLS variants of each.
#'
#' @param x Predictor data: an \eqn{n \times p}{n x p} numeric matrix,
#'   data.frame, tibble, or tsibble.
#'   For data.frame/tibble/tsibble inputs, \code{y} can be a column name
#'   and date columns are auto-detected.
#' @param y Response vector (length \eqn{n}), or a column name string
#'   when \code{x} is a data.frame/tibble/tsibble.
#' @param method Estimation method. One of:
#'   \describe{
#'     \item{\code{"lasso"}}{Standard Lasso (Tibshirani, 1996).
#'       \eqn{\ell_1}{L1}-penalized least squares.}
#'     \item{\code{"alasso"}}{Adaptive Lasso (Zou, 2006).
#'       Uses data-dependent penalty weights
#'       \eqn{w_j = 1/|\hat\beta_j^{init}|^\gamma}{w_j = 1/|beta_init_j|^gamma}
#'       from an initial estimator.}
#'     \item{\code{"talasso"}}{Twin Adaptive Lasso (Lee, Shi, and Gao, 2022).
#'       Two-round adaptive Lasso: fits ALasso first, then uses its coefficients
#'       as weights in a second penalized regression.}
#'     \item{\code{"post_lasso"}, \code{"post_alasso"}, \code{"post_talasso"}}{
#'       Post-selection OLS: runs the corresponding Lasso method to select
#'       variables, then re-estimates by OLS on the selected set to reduce
#'       shrinkage bias.}
#'   }
#' @param predictive Logical. If \code{TRUE}, automatically lag-align the data
#'   so that row \eqn{t} pairs \eqn{y_t} with \eqn{x_{t-1}}. If \code{FALSE},
#'   the data is used as-is. Default \code{NULL} auto-detects: \code{TRUE} when
#'   date information is present, \code{FALSE} otherwise.
#' @param date_col Optional date column name (for data.frame/tibble input).
#' @param ... Additional arguments forwarded to \code{\link{train_lasso}}
#'   (tuning) and the underlying fitting functions.
#'   Key arguments:
#'   \describe{
#'     \item{\code{train_method}}{Tuning parameter selection method
#'       (default \code{"timeslice"}).
#'       Options: \code{"timeslice"}, \code{"cv"}, \code{"cv_random"},
#'       \code{"aic"}, \code{"bic"}, \code{"aicc"}, \code{"hqc"}.}
#'     \item{\code{scale_x}}{Logical; standardize predictors before fitting.
#'       Default \code{TRUE}, except for \code{method = "alasso"} when the
#'       initial estimator resolves to OLS (low-dimensional case), where it
#'       defaults to \code{FALSE} so the weighted Lasso runs on the raw scale
#'       as in Lee, Shi, and Gao (2022). TALasso inherits this in its first
#'       round; its second round defaults to \code{FALSE}.}
#'     \item{\code{intercept}}{Logical; include an intercept
#'       (default \code{TRUE}).}
#'     \item{\code{gamma}}{Exponent for adaptive weights in ALasso/TALasso
#'       (default 1).}
#'     \item{\code{model}}{Initial estimator for adaptive weights:
#'       \code{"ols"}, \code{"lasso"}, \code{"ridge"}, or \code{NULL}
#'       (auto-selects based on \eqn{n} vs. \eqn{p}; default \code{NULL}).}
#'     \item{\code{lambda_init}}{Penalty for the Lasso/ridge \emph{initial}
#'       estimator generating the adaptive weights (ALasso, and TALasso's
#'       first round). Default \code{NULL}: tuned internally by 10-fold block
#'       CV. Ignored when the initial estimator is OLS. The same value is
#'       used in tuning and estimation, so the two stay coherent.}
#'     \item{\code{nlambda}}{Number of candidate \eqn{\lambda} values
#'       (default 100).}
#'     \item{\code{k}}{Number of folds for CV methods (default 10).}
#'   }
#'   See \code{\link{train_lasso}} for the full list of tuning parameters.
#'
#' @return An object of class \code{"lasforecast_model"} (a list) with:
#'   \describe{
#'     \item{\code{fit}}{The fitted model object (\code{glmnet} for Lasso,
#'       internal list for ALasso/TALasso, \code{NULL} for Post-Lasso).}
#'     \item{\code{lambda}}{Selected penalty parameter.}
#'     \item{\code{coefficients}}{Numeric vector of estimated coefficients
#'       (intercept first, then slopes).}
#'     \item{\code{method}}{Character string identifying the method used.}
#'   }
#'   S3 methods: \code{\link{predict.lasforecast_model}},
#'   \code{\link[=print.lasforecast_model]{print}},
#'   \code{\link[=coef.lasforecast_model]{coef}}.
#'
#' @references
#' Lee, J. H., Shi, Z., and Gao, Z. (2022).
#' "On LASSO for Predictive Regression."
#' \emph{Journal of Econometrics}, 229(2), 322--349.
#'
#' Mei, Z., & Shi, Z. (2024).
#' "On LASSO for High Dimensional Predictive Regression."
#' \emph{Journal of Econometrics}, 242(2): 105809.
#'
#' @seealso \code{\link{roll_predict}} for rolling-window forecasting,
#'   \code{\link{train_lasso}} for tuning parameter details,
#'   \code{\link{xdlasso}} for debiased inference.
#'
#' @examples
#' x <- matrix(rnorm(100 * 5), 100, 5)
#' y <- x[, 1] + rnorm(100)
#'
#' # Lasso with BIC tuning
#' mod <- lasso(x, y, method = "lasso", train_method = "bic")
#' coef(mod)
#'
#' # Adaptive Lasso
#' mod_a <- lasso(x, y, method = "alasso")
#'
#' # Post-Lasso OLS (Lasso selection, then OLS on selected set)
#' mod_pl <- lasso(x, y, method = "post_lasso")
#'
#' @export
lasso <- function(x, y = NULL,
                  method = c(
                      "lasso", "alasso", "talasso",
                      "post_lasso", "post_alasso", "post_talasso"
                  ),
                  predictive = NULL, date_col = NULL, ...) {
    method <- match.arg(method)
    inp <- prepare_input(x, y, date_col)
    x <- inp$x_mat
    y <- inp$y_vec

    if (is.null(predictive)) predictive <- !is.null(inp$dates)
    if (predictive) {
        n_raw <- nrow(x)
        x <- x[-n_raw, , drop = FALSE]
        y <- y[-1]
    }
    args <- list(...)

    if (startsWith(method, "post_")) {
        base_method <- sub("^post_", "", method)
        est <- lasso(x, y, method = base_method, ...)
        coef_post <- post_lasso(x, y, est$coefficients)
        return(structure(
            list(fit = NULL, coefficients = coef_post, method = method),
            class = "lasforecast_model"
        ))
    }

    intercept <- args$intercept %||% TRUE
    scale_x <- args$scale_x %||% TRUE
    gamma <- args$gamma %||% 1

    # Default scale convention follows the initial estimator (Lee, Shi & Gao,
    # 2022): with OLS weights the penalty acts on the raw scale (fixed-p
    # theory); with a shrinkage initial estimator (high-dimensional case) the
    # standardized convention applies. An explicit scale_x always wins.
    if (method == "alasso" && is.null(args$scale_x)) {
        init_model <- args$model %||% (if (ncol(x) > 2 * sqrt(nrow(x))) "lasso" else "ols")
        scale_x <- (init_model != "ols")
    }

    if (method %in% c("lasso", "alasso")) {
        ada <- (method == "alasso")

        train_args <- args[names(args) %in% names(formals(train_lasso))]
        train_args$x <- x
        train_args$y <- y
        train_args$ada <- ada
        train_args$intercept <- intercept
        train_args$scale_x <- scale_x
        if (ada) train_args$gamma <- gamma
        lambda_opt <- do.call(train_lasso, train_args)

        if (ada) {
            fit <- adalasso(x, y,
                lambda = lambda_opt, lambda_init = args$lambda_init,
                intercept = intercept,
                scale_x = scale_x, gamma = gamma, model = args$model
            )
            coefs <- c(fit$ahat, fit$bhat)
        } else {
            glmnet_args <- args[names(args) %in% names(formals(glmnet::glmnet))]
            glmnet_args$x <- x
            glmnet_args$y <- y
            glmnet_args$lambda <- lambda_opt
            glmnet_args$intercept <- intercept
            glmnet_args$standardize <- scale_x
            fit <- do.call(glmnet::glmnet, glmnet_args)
            coefs <- as.numeric(coef(fit))
        }

        structure(
            list(fit = fit, lambda = lambda_opt, coefficients = coefs, method = method),
            class = "lasforecast_model"
        )
    } else if (method == "talasso") {
        est_ada <- lasso(x, y, method = "alasso", ...)
        if (!identical(est_ada$fit$model, "ols")) {
            warning("TALasso is designed for low-dimensional predictive regressions, but the first-stage adaptive Lasso used a shrinkage initial estimator ('", est_ada$fit$model, "', chosen when p > 2 * sqrt(n) or set via model). Consider method = 'alasso' for high-dimensional data.")
        }
        b_first <- est_ada$coefficients[-1]

        train_args <- args[names(args) %in% names(formals(train_talasso))]
        train_args$x <- x
        train_args$y <- y
        train_args$b_first <- b_first
        if (is.null(train_args$gamma)) train_args$gamma <- gamma
        if (is.null(train_args$intercept)) train_args$intercept <- intercept
        if (is.null(train_args$scale_x)) train_args$scale_x <- FALSE
        lambda_ta <- do.call(train_talasso, train_args)

        fit <- talasso(x, y, b_first,
            lambda = lambda_ta,
            gamma = train_args$gamma, intercept = train_args$intercept,
            scale_x = train_args$scale_x
        )

        structure(
            list(
                fit = fit, lambda = lambda_ta,
                coefficients = c(fit$ahat, fit$bhat), method = "talasso"
            ),
            class = "lasforecast_model"
        )
    }
}


#' IVX-Desparsified Lasso Inference for Predictive Regression
#'
#' Performs inference on focal predictors in a high-dimensional predictive
#' regression using the IVX-desparsified Lasso (XDlasso) method of
#' Gao, Lee, Mei, and Shi (2026). Provides asymptotically valid \eqn{t}-tests
#' and confidence intervals that are robust to both Lasso regularization bias
#' and Stambaugh bias from persistent regressors.
#'
#' The model is:
#' \deqn{y_t = d_{t-1}'\theta + x_{t-1}'\beta + u_t}
#' where \eqn{d_t} are the focal predictors of interest and \eqn{x_t} are
#' high-dimensional controls selected by Lasso. The method constructs IVX-type
#' instruments to debias the Lasso estimate of \eqn{\theta}, yielding standard
#' normal inference.
#'
#' @param w Predictor data. Can be an \eqn{n \times p}{n x p} numeric matrix,
#'   data.frame, tibble, or tsibble containing \emph{all} regressors (focal
#'   predictors and controls). For data.frame/tibble/tsibble inputs, \code{y}
#'   can be a column name and date columns are auto-detected.
#' @param y Response vector (length \eqn{n}), or a column name string
#'   when \code{w} is a data.frame/tibble/tsibble.
#' @param d Character vector of focal variable names (column names in
#'   \code{w}). Alternative to \code{d_ind}; exactly one must be specified.
#' @param d_ind Integer vector of column indices identifying the focal
#'   predictors. Alternative to \code{d}.
#' @param predictive Logical. If \code{TRUE}, automatically lag-align the data
#'   so that row \eqn{t} pairs \eqn{y_t} with \eqn{w_{t-1}}. If \code{FALSE},
#'   the data is used as-is (the user is responsible for alignment). Default
#'   \code{NULL} auto-detects: \code{TRUE} when date information is present,
#'   \code{FALSE} otherwise.
#' @param date_col Optional date column name (for data.frame/tibble input).
#' @param intercept Logical. Include an intercept term (default \code{TRUE}).
#' @param ... Arguments passed to \code{\link{debias_ivx}}.
#'   Key arguments:
#'   \describe{
#'     \item{\code{standardize}}{Logical; standardize regressors before
#'       Lasso fitting (default \code{TRUE}).}
#'     \item{\code{c_z}, \code{a}}{IVX instrument tuning parameters
#'       (defaults 5 and 0.5).}
#'     \item{\code{train_method}}{Lambda selection for internal Lasso fits
#'       (default \code{"timeslice"}). Same options as \code{\link{train_lasso}}.}
#'     \item{\code{se_type}}{Standard error type: \code{"iid"} (default),
#'       \code{"robust"}, or \code{"HAC"}.}
#'     \item{\code{joint_test}}{Logical; compute a joint Wald test for all
#'       focal coefficients (default \code{FALSE}).}
#'   }
#'   See \code{\link{debias_ivx}} for the full parameter list.
#'
#' @return An object of class \code{"lasforecast_debias_ivx"} (a list) with:
#'   \describe{
#'     \item{\code{theta_hat_las}}{Lasso estimates of the focal coefficients.}
#'     \item{\code{theta_hat_ivx}}{Debiased IVX estimates of \eqn{\theta}
#'       (the main inferential quantities).}
#'     \item{\code{sigma_hat_ivx}}{Standard errors for \code{theta_hat_ivx}.}
#'     \item{\code{d_ind}}{Focal predictor indices.}
#'     \item{\code{term_names}}{Predictor names (from \code{colnames(w)}).}
#'   }
#'   S3 methods: \code{\link[=print.lasforecast_debias_ivx]{print}},
#'   \code{\link[=summary.lasforecast_debias_ivx]{summary}} (coefficient
#'   table with \eqn{t}-statistics and \eqn{p}-values),
#'   \code{\link[=coef.lasforecast_debias_ivx]{coef}}.
#'
#' @references
#' Gao, Z., Lee, J. H., Mei, Z., and Shi, Z. (2026).
#' "LASSO Inference for High Dimensional Predictive Regressions."
#'
#' @seealso \code{\link{lasso}} for point estimation and forecasting,
#'   \code{\link{debias_ivx}} for the underlying estimation engine.
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 10
#' w <- matrix(rnorm(n * p), n, p)
#' colnames(w) <- paste0("x", 1:p)
#' y <- w[, 1] * 0.5 + rnorm(n)
#'
#' fit <- xdlasso(w, y, d_ind = 1)
#' summary(fit)
#'
#' @export
xdlasso <- function(w, y = NULL, d = NULL, d_ind = NULL,
                    predictive = NULL, date_col = NULL, intercept = TRUE, ...) {
    inp <- prepare_input(w, y, date_col)
    w_mat <- inp$x_mat
    y_vec <- inp$y_vec
    dates <- inp$dates

    if (!is.null(d) && !is.null(d_ind)) {
        stop("Specify either 'd' (names) or 'd_ind' (indices), not both.")
    }
    if (is.null(d) && is.null(d_ind)) {
        stop("One of 'd' (focal variable names) or 'd_ind' (focal variable indices) must be specified.")
    }

    if (!is.null(d)) {
        d_ind <- match(d, colnames(w_mat))
        if (anyNA(d_ind)) {
            stop("Focal variable(s) not found: ", paste(d[is.na(d_ind)], collapse = ", "))
        }
    }

    if (is.null(predictive)) predictive <- !is.null(dates)
    if (predictive) {
        n_raw <- nrow(w_mat)
        w_mat <- w_mat[-n_raw, , drop = FALSE]
        y_vec <- y_vec[-1]
    }

    fit <- debias_ivx(w = w_mat, y = y_vec, d_ind = d_ind, intercept = intercept, ...)

    term_names <- colnames(w_mat)
    if (is.null(term_names)) term_names <- paste0("V", seq_len(ncol(w_mat)))

    fit$d_ind <- d_ind
    fit$term_names <- term_names

    class(fit) <- "lasforecast_debias_ivx"
    fit
}


#' Estimate OLS Model
#'
#' Fits an OLS model.
#'
#' @param x Predictor data. Can be a matrix, data.frame, tibble, or tsibble.
#' @param y Response vector, or column name string when \code{x} is a data frame.
#' @param intercept Logical
#' @param date_col Optional date column name (for data.frame/tibble input).
#' @param ... Additional arguments (ignored).
#'
#' @return A list with class "lasforecast_model".
#' @export
ols <- function(x, y = NULL, intercept = TRUE, date_col = NULL, ...) {
    inp <- prepare_input(x, y, date_col)
    x <- inp$x_mat
    y <- inp$y_vec
    fit <- stats::lsfit(x, y, intercept = intercept)

    structure(
        list(
            fit = fit,
            coefficients = as.numeric(fit$coefficients),
            method = "ols"
        ),
        class = "lasforecast_model"
    )
}

#' Estimate Random Walk
#' @param x Predictor data. Can be a matrix, data.frame, tibble, or tsibble.
#' @param y Response vector, or column name string when \code{x} is a data frame.
#' @param date_col Optional date column name.
#' @param ... Additional arguments (ignored).
#' @param type Type of forecast target y. This is only relevant for RW method.
#'  "diff" for returns/differences and RW predicts zeros; "level" for levels and RW predicts the previous value. Default is "level".
#' @export
rw <- function(x, y = NULL, date_col = NULL, type = "level", ...) {
    inp <- prepare_input(x, y, date_col)
    x <- inp$x_mat
    y <- inp$y_vec

    if (type == "diff") {
        coef_rw <- rep(0, ncol(x) + 1)
    } else if (type == "level") {
        coef_rw <- c(tail(y, 1), rep(0, ncol(x)))
    } else {
        stop("Invalid type. Choose 'diff' or 'level'.")
    }
    structure(
        list(
            type = type,
            coefficients = coef_rw,
            method = "rw"
        ),
        class = "lasforecast_model"
    )
}

#' Estimate RW with Drift
#' @param x Predictor data. Can be a matrix, data.frame, tibble, or tsibble.
#' @param y Response vector, or column name string when \code{x} is a data frame.
#' @param date_col Optional date column name.
#' @param ... Additional arguments (ignored).
#' @export
rwwd <- function(x, y = NULL, date_col = NULL, ...) {
    inp <- prepare_input(x, y, date_col)
    x <- inp$x_mat
    y <- inp$y_vec
    mu <- mean(y)
    structure(
        list(
            coefficients = c(mu, rep(0, ncol(x))),
            method = "rwwd"
        ),
        class = "lasforecast_model"
    )
}

#' Estimate Univariate ARMA Benchmark
#'
#' Fits a univariate ARMA(\eqn{p}, \eqn{d}, \eqn{q}) model on the response
#' \code{y} and returns an \eqn{h}-step-ahead forecast. This is a time-series
#' benchmark, parallel to \code{\link{rw}} and \code{\link{rwwd}}: the
#' predictors \code{x} are \emph{not} used. When \code{order} is \code{NULL},
#' the AR and MA orders are selected automatically by information criterion
#' over a grid (the differencing order \code{d} is held fixed).
#'
#' @param x Optional predictor data, accepted only for interface symmetry with
#'   the other benchmark estimators; it is ignored. When \code{x} is a
#'   data.frame/tibble/tsibble, \code{y} may be given as a column name.
#' @param y Response vector (length \eqn{n}), or a column name string when
#'   \code{x} is a data.frame/tibble/tsibble.
#' @param order Optional ARMA order \code{c(p, d, q)}. If \code{NULL} (default),
#'   \eqn{p} and \eqn{q} are chosen automatically by \code{ic} over the grid
#'   \eqn{\{0, \dots, }\code{max_p}\eqn{\} \times \{0, \dots, }\code{max_q}\eqn{\}}.
#' @param h Forecast horizon (default 1). The function returns the
#'   \eqn{h}-step-ahead forecast from the end of \code{y}.
#' @param include.mean Logical; include a mean term (default \code{TRUE}).
#'   Ignored by \code{\link[stats]{arima}} when \code{d > 0}.
#' @param d Differencing order used during automatic selection (default 0, i.e.
#'   treat \code{y} as stationary). Ignored when \code{order} is supplied.
#' @param max_p,max_q Maximum AR and MA orders searched during automatic
#'   selection (defaults 5).
#' @param ic Information criterion for automatic selection: \code{"aic"}
#'   (default) or \code{"bic"}.
#' @param arma_method Fitting method passed to \code{\link[stats]{arima}}:
#'   \code{"CSS-ML"} (default), \code{"ML"}, or \code{"CSS"}.
#' @param date_col Optional date column name (for data.frame/tibble input).
#' @param ... Additional arguments (ignored).
#'
#' @return An object of class \code{"lasforecast_model"} with components
#'   \code{fit} (the \code{\link[stats]{arima}} object), \code{forecast} (the
#'   \eqn{h}-step point forecast), \code{arma_order}, \code{df} (number of
#'   estimated parameters), \code{coefficients} (\code{NA}; ARMA has no
#'   predictor coefficients), and \code{method = "arma"}.
#'
#' @seealso \code{\link{rw}}, \code{\link{rwwd}}, \code{\link{ols}},
#'   \code{\link{roll_predict}}, \code{\link{backtest}}.
#'
#' @examples
#' set.seed(1)
#' y <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.3), n = 120))
#' arma(y = y, order = c(1, 0, 1))   # fixed order
#' arma(y = y)                        # automatic order selection
#'
#' @importFrom stats arima
#' @export
arma <- function(x = NULL, y = NULL, order = NULL, h = 1,
                 include.mean = TRUE, d = 0,
                 max_p = 5, max_q = 5, ic = c("aic", "bic"),
                 arma_method = "CSS-ML", date_col = NULL, ...) {
    ic <- match.arg(ic)

    # Univariate: extract y only. x is accepted for interface symmetry but ignored.
    if (is.null(x)) {
        y <- as.numeric(y)
    } else {
        y <- prepare_input(x, y, date_col)$y_vec
    }

    # Fitting closure: returns an arima fit or NULL on failure (rolling windows
    # occasionally hit non-convergence; warnings are suppressed to avoid flooding).
    fit_arima <- function(ord, mth = arma_method) {
        tryCatch(
            suppressWarnings(stats::arima(y, order = ord,
                                          include.mean = include.mean, method = mth)),
            error = function(e) NULL
        )
    }

    if (is.null(order)) {
        sel <- arma_auto_select(y, fit_arima, d = d, max_p = max_p,
                                max_q = max_q, ic = ic)
        fit <- sel$fit
        order <- sel$order
    } else {
        fit <- fit_arima(order)
    }

    # Robust fallbacks so a benchmark never returns a corrupting forecast.
    if (is.null(fit)) {
        fit <- fit_arima(c(0, d, 0), mth = "ML")
        order <- c(0, d, 0)
    }
    if (is.null(fit)) {
        mu <- mean(y)
        return(structure(
            list(fit = NULL, forecast = mu, arma_order = c(0, 0, 0),
                 df = 1, coefficients = NA_real_, method = "arma"),
            class = "lasforecast_model"
        ))
    }

    fc <- as.numeric(predict(fit, n.ahead = h)$pred)[h]

    structure(
        list(
            fit = fit,
            forecast = fc,
            arma_order = order,
            df = length(fit$coef),
            coefficients = NA_real_,
            method = "arma"
        ),
        class = "lasforecast_model"
    )
}

#' Automatic ARMA order selection by information criterion
#'
#' Internal grid search over AR and MA orders (with \code{d} held fixed) used by
#' \code{\link{arma}} when \code{order} is not supplied. The criterion is built
#' from the maximized log-likelihood for comparability across orders:
#' \eqn{-2\ell + 2k} (AIC) or \eqn{-2\ell + \log(n)\,k} (BIC), where \eqn{k} is
#' the number of estimated parameters including the innovation variance.
#'
#' @param y Numeric response vector.
#' @param fit_arima A fitting closure returning an \code{arima} object or
#'   \code{NULL} on failure.
#' @param d Differencing order (held fixed).
#' @param max_p,max_q Maximum AR and MA orders.
#' @param ic \code{"aic"} or \code{"bic"}.
#'
#' @return A list with \code{fit} (best \code{arima} object, or \code{NULL} if
#'   every candidate failed) and \code{order} (the selected \code{c(p, d, q)}).
#'
#' @keywords internal
arma_auto_select <- function(y, fit_arima, d = 0, max_p = 5, max_q = 5, ic = "aic") {
    best_ic <- Inf
    best_fit <- NULL
    best_order <- c(0, d, 0)
    for (p in 0:max_p) {
        for (q in 0:max_q) {
            fit <- fit_arima(c(p, d, q))
            if (is.null(fit) || !is.finite(fit$loglik)) next
            k <- length(fit$coef) + 1 # + innovation variance
            val <- if (ic == "aic") {
                -2 * fit$loglik + 2 * k
            } else {
                -2 * fit$loglik + log(fit$nobs) * k
            }
            if (is.finite(val) && val < best_ic) {
                best_ic <- val
                best_fit <- fit
                best_order <- c(p, d, q)
            }
        }
    }
    list(fit = best_fit, order = best_order)
}

#' Predict method for LasForecast models
#'
#' @param object Model object
#' @param newx New data matrix
#' @param ... Additional arguments
#'
#' @export
predict.lasforecast_model <- function(object, newx, ...) {
    # ARMA is univariate: the h-step forecast is precomputed at fit time and
    # does not depend on newx.
    if (identical(object$method, "arma")) {
        return(object$forecast)
    }

    if (is.vector(newx)) {
        newx <- matrix(newx, nrow = 1)
    }

    coefs <- object$coefficients

    if (length(coefs) == ncol(newx) + 1) {
        y_hat <- as.vector(cbind(1, newx) %*% coefs)
    } else if (length(coefs) == ncol(newx)) {
        y_hat <- as.vector(newx %*% coefs)
    } else {
        stop(
            "Dimension mismatch: Coefficients length ", length(coefs),
            " vs New Data columns ", ncol(newx)
        )
    }

    y_hat
}
