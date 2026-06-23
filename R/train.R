# ----
#' Tuning parameter selection for LASSO
#'
#' Selects the penalty parameter lambda for LASSO or Adaptive LASSO
#' via cross-validation, time-series CV, or information criteria.
#'
#' @param x Predictor matrix (\eqn{n \times p}).
#' @param y Response vector (length \eqn{n}).
#' @param ada Logical. If \code{TRUE} (default), tune for Adaptive LASSO;
#'   if \code{FALSE}, tune for standard LASSO.
#' @param gamma Exponent for adaptive penalty weights
#'   \eqn{\hat\tau_j = |\hat\theta_j^{init}|^{-\gamma}}. Default 1.
#' @param alpha Elastic-net mixing parameter passed to \code{glmnet}.
#'   Use 1 for LASSO (default) and 0 for ridge.
#' @param intercept Logical. Include an intercept term (default \code{TRUE}).
#' @param model Initial estimator for adaptive weights: \code{"lasso"},
#'   \code{"ridge"}, or \code{"ols"}. Default \code{NULL} auto-selects
#'   based on dimensionality (\eqn{n} vs \eqn{p}).
#' @param lambda_init Penalty for the Lasso/ridge \emph{initial} estimator
#'   generating the adaptive weights when \code{ada = TRUE}. Default
#'   \code{NULL}: tuned internally by 10-fold block CV. Ignored when
#'   \code{ada = FALSE} or the initial estimator is OLS.
#' @param scale_x Logical. Standardize predictors before fitting (default \code{FALSE}).
#' @param lambda_seq Candidate sequence of lambda values. If \code{NULL}
#'   (default), the function generates the sequence automatically.
#' @param train_method Tuning parameter selection method. One of:
#'   \describe{
#'     \item{\code{"cv"}}{Block cross-validation. Splits the sample into \code{k}
#'       contiguous blocks; each block serves as the validation fold in turn.
#'       \preformatted{
#'  Iter 1: train train train train  val
#'  Iter 2: train train train  val  train
#'  Iter 3: train train  val  train train
#'       }}
#'     \item{\code{"cv_random"}}{Random cross-validation. Assigns observations to
#'       \code{k} folds randomly, ignoring temporal order. Appropriate only when
#'       serial dependence is weak.}
#'     \item{\code{"timeslice"}}{Time-series cross-validation. Training always uses
#'       past data; validation uses future data. Controlled by \code{initial_window},
#'       \code{horizon}, \code{fixed_window}, and \code{skip}.
#'       \preformatted{
#'  Rolling (fixed_window = TRUE):
#'  Iter 1: ===train===|=val=|
#'  Iter 2:   ===train===|=val=|
#'  Iter 3:     ===train===|=val=|
#'
#'  Expanding (fixed_window = FALSE):
#'  Iter 1: ===train===|=val=|
#'  Iter 2: ====train====|=val=|
#'  Iter 3: =====train=====|=val=|
#'       }}
#'     \item{\code{"aic"}, \code{"bic"}, \code{"aicc"}, \code{"hqc"}}{Information
#'       criterion. Fits the model along the full lambda path and selects the lambda
#'       minimizing the chosen criterion. Fast; no sample splitting required.}
#'   }
#' @param nlambda Number of candidate lambda values (default 100).
#' @param lambda_min_ratio Ratio determining the smallest lambda:
#'   \code{lambda_min = lambda_min_ratio * lambda_max} (default 0.0001).
#' @param k Number of folds for \code{"cv"} and \code{"cv_random"} methods
#'   (default 10).
#' @param initial_window Size of the first training set for \code{"timeslice"}
#'   (default: 70\% of \eqn{n}).
#' @param horizon Number of periods in each validation set for
#'   \code{"timeslice"} (default 1).
#' @param fixed_window Logical. If \code{TRUE} (default), the training window
#'   rolls forward at fixed size; if \code{FALSE}, it expands.
#' @param skip Number of periods to skip between successive splits in
#'   \code{"timeslice"} (default 0).
#'
#' @return The selected lambda value (numeric scalar).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(50 * 5), 50, 5)
#' y <- rnorm(50)
#' train_lasso(x, y)
#' }
train_lasso <- function(
  x,
  y,
  ada = TRUE,
  gamma = 1,
  model = NULL,
  lambda_init = NULL,
  alpha = 1,
  intercept = TRUE,
  scale_x = FALSE,
  lambda_seq = NULL,
  train_method = "timeslice",
  nlambda = 100,
  lambda_min_ratio = 0.0001,
  k = 10,
  initial_window = ceiling(nrow(x) * 0.7),
  horizon = 1,
  fixed_window = TRUE,
  skip = 0
) {
    n <- nrow(x)
    p <- ncol(x)

    if (ada && alpha != 1) {
        stop("Adaptive lasso tuning only supports alpha = 1. Use ada = FALSE for ridge tuning.")
    }

    w <- NULL
    if (ada) {
        init_result <- init_est(x, y, lambda_init = lambda_init, gamma = gamma,
                                intercept = intercept, scale_x = scale_x, model = model)
        w <- init_result$w
        coef_init <- init_result$coef_init
    }

    if (is.null(lambda_seq)) {
        if (alpha == 1) {
            if (ada) {
                # Coordinate j enters the all-zero solution at
                # lambda = |corr_j| / w_j = |corr_j| * |theta_init_j|^gamma,
                # so the grid top is the max over j of these per-coordinate
                # products. max-of-products (not the old product-of-maxes upper
                # bound) keeps lambda_max equal to the exact entry threshold and
                # scale-equivariant at gamma = 1 (Lee, Shi & Gao 2022). corr_j
                # is read on the same operating scale as the weights.
                cross_j <- get_lasso_lambda_max(x, y, scale_x = scale_x,
                                                intercept = intercept, per_coord = TRUE)
                coef_slope <- if (intercept) coef_init[-1] else coef_init
                if (scale_x) coef_slope <- coef_slope * apply(x, 2, sd_n)
                lambda_max <- max(cross_j * abs(coef_slope)^gamma)
                # Degenerate case: the initial estimator screened out every
                # variable, so all weights sit at the floor and any lambda
                # yields the empty model; fall back to the plain-Lasso cap so
                # the log-spaced grid stays defined.
                if (lambda_max == 0) lambda_max <- max(cross_j)
            } else {
                lambda_max <- get_lasso_lambda_max(x, y, scale_x = scale_x,
                                                   intercept = intercept)
            }
            # Nudge the grid top just above the exact entry threshold so the
            # intercept-only model is reliably the first grid point: glmnet's
            # reported df can flip 0 <-> 1 exactly at the boundary under
            # predictor rescaling, which would otherwise drop the empty model
            # from the path and bias IC selection (breaking the scale
            # invariance of standardized Lasso).
            lambda_max <- lambda_max * (1 + 1e-7)
            lambda_seq <- get_lambda_seq(lambda_max, lambda_min_ratio = lambda_min_ratio, nlambda = nlambda)
        } else {
            lambda_seq <- glmnet::glmnet(x = x, y = y, alpha = alpha, intercept = intercept,
                                         standardize = scale_x, nlambda = nlambda,
                                         lambda.min.ratio = lambda_min_ratio)$lambda
        }
    }

    glmnet_args <- list(x = x, y = y, lambda = lambda_seq, alpha = alpha,
                        intercept = intercept, standardize = scale_x, nfolds = k)
    if (ada) {
        glmnet_args$penalty.factor <- w
        glmnet_args$lambda <- lambda_seq * sum(w) / p
    }

    if (train_method %in% c("cv", "cv_random")) {
        .train_cv(glmnet_args, train_method, n, k, ada, p, w)
    } else if (train_method %in% c("aic", "bic", "aicc", "hqc")) {
        .train_ic(glmnet_args, train_method, lambda_seq, n, intercept, alpha, ada, p, w)
    } else if (train_method == "timeslice") {
        .train_timeslice(x, y, alpha, intercept, scale_x, lambda_seq,
                         initial_window, horizon, fixed_window, skip, ada, w, p)
    } else {
        stop("Invalid train_method input.")
    }
}


.train_cv <- function(glmnet_args, train_method, n, k, ada, p, w) {
    if (train_method == "cv") {
        glmnet_args$foldid <- foldid_vec(n, k = k)
    }
    cv_las <- do.call(glmnet::cv.glmnet, glmnet_args)
    lambda_cv <- cv_las$lambda.min
    if (ada) lambda_cv * p / sum(w) else lambda_cv
}


.train_ic <- function(glmnet_args, train_method, lambda_seq, n, intercept, alpha, ada, p, w) {
    if (alpha == 0) {
        stop("Information-criterion tuning is not supported for ridge. Use cv, cv_random, or timeslice instead.")
    }
    glm_est <- do.call(glmnet::glmnet, glmnet_args)
    lambda_fit <- glm_est$lambda

    y_hat <- as.matrix(predict(glm_est, newx = glmnet_args$x))
    e_hat <- matrix(glmnet_args$y, n, length(lambda_fit)) - y_hat
    mse <- colMeans(e_hat^2)

    nvar <- glm_est$df + intercept
    ic <- switch(train_method,
        "aic" = n * log(mse) + 2 * nvar,
        "bic" = n * log(mse) + log(n) * nvar,
        "aicc" = n * log(mse) + 2 * nvar + 2 * nvar * (nvar + 1) / (n - nvar - 1),
        "hqc" = n * log(mse) + 2 * nvar * log(log(n))
    )

    # Guard the saturated tail: mse = 0 gives log(mse) = -Inf and the AICc
    # correction divides by (n - nvar - 1) <= 0, either of which would make
    # which.min lock onto a degenerate (over-fit) lambda. Exclude them.
    ic[!is.finite(ic)] <- Inf
    lambda_temp <- lambda_fit[which.min(ic)]
    if (ada) lambda_temp * p / sum(w) else lambda_temp
}


.train_timeslice <- function(x, y, alpha, intercept, scale_x, lambda_seq,
                              initial_window, horizon, fixed_window, skip, ada, w, p) {
    n <- nrow(x)
    slices <- create_time_slices(n, initial_window, horizon, fixed_window, skip)

    fit_lambda <- if (ada) lambda_seq * sum(w) / p else lambda_seq

    glmnet_base <- list(alpha = alpha, intercept = intercept, standardize = scale_x)
    if (ada) glmnet_base$penalty.factor <- w

    mse_sum <- rep(0, length(fit_lambda))
    for (j in seq_along(slices$train)) {
        args <- c(list(x = x[slices$train[[j]], , drop = FALSE],
                       y = y[slices$train[[j]]],
                       lambda = fit_lambda), glmnet_base)
        fit <- do.call(glmnet::glmnet, args)
        pred <- predict(fit, newx = x[slices$test[[j]], , drop = FALSE])
        mse_sum <- mse_sum + colMeans((y[slices$test[[j]]] - pred)^2)
    }

    best_lambda <- fit_lambda[which.min(mse_sum)]
    if (ada) best_lambda * p / sum(w) else best_lambda
}


create_time_slices <- function(n, initial_window, horizon, fixed_window = TRUE, skip = 0) {
    stops <- seq(initial_window, n - horizon, by = skip + 1)
    train_list <- vector("list", length(stops))
    test_list <- vector("list", length(stops))
    for (i in seq_along(stops)) {
        start <- if (fixed_window) stops[i] - initial_window + 1 else 1
        train_list[[i]] <- start:stops[i]
        test_list[[i]] <- (stops[i] + 1):(stops[i] + horizon)
    }
    list(train = train_list, test = test_list)
}


foldid_vec <- function(TT, k) {
    seq_interval <- split(1:TT, ceiling(seq_along(1:TT) / (TT / k)))
    id <- rep(0, TT)
    for (j in 1:k) {
        id[seq_interval[[j]]] <- j
    }
    id
}

# -----

#' Do parameter tuning for Twin Adaptive Lasso (TALasso)
#'
#' If all variables are killed in the first step: return NA\cr
#' If more than 1 variables are left: just repeat the training process for alasso\cr
#' If only 1 variable remained: use a brute-force process do the cross-validation.\cr
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param b_first estimated slope from first step alasso
#' @param gamma Parameter controlling the inverse of first step estimate. By default = 1.
#' @param model Ignored. The second-step initial estimator is always OLS,
#'   because the first-step-selected set is low-dimensional by construction and
#'   talasso() builds its weights from the OLS fit. Retained for interface
#'   compatibility with \code{train_lasso}.
#' @param intercept A boolean: include an intercept term or not
#' @param scale_x A boolean: standardize the design matrix or not
#' @param lambda_seq Candidate sequnece of parameters. If NULL, the function generates the sequnce.
#' @param train_method "timeslice", "cv",  "aic", "bic", "aicc", "hqc"
#' @param nlambda Number of candidate lambda values (default 100).
#' @param lambda_min_ratio Ratio determining the smallest lambda:
#'   \code{lambda_min = lambda_min_ratio * lambda_max} (default 0.0001).
#' @param k Number of folds for \code{"cv"} and \code{"cv_random"} methods
#'   (default 10).
#' @param initial_window Size of the first training set for \code{"timeslice"}
#'   (default: 70\% of \eqn{n}).
#' @param horizon Number of periods in each validation set for
#'   \code{"timeslice"} (default 1).
#' @param fixed_window Logical. If \code{TRUE} (default), the training window
#'   rolls forward at fixed size; if \code{FALSE}, it expands.
#' @param skip Number of periods to skip between successive splits in
#'   \code{"timeslice"} (default 0).
#'
#' @return The selected lambda value (numeric scalar).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(50 * 5), 50, 5)
#' y <- rnorm(50)
#' b_first <- rnorm(5)
#' train_talasso(x, y, b_first = b_first)
#' }
train_talasso <- function(x,
                          y,
                          b_first,
                          gamma = 1,
                          model = NULL,
                          intercept = TRUE,
                          scale_x = FALSE,
                          train_method = "timeslice",
                          lambda_seq = NULL,
                          nlambda = 100,
                          lambda_min_ratio = 0.0001,
                          k = 10,
                          initial_window = ceiling(nrow(x) * 0.7),
                          horizon = 1,
                          fixed_window = TRUE,
                          skip = 0) {
    n <- nrow(x)
    p <- ncol(x)

    selected <- (b_first != 0)
    p_selected <- sum(selected)
    xx <- as.matrix(x[, selected])

    if (p_selected == 0) {
        warning("All predictors screened out in first step. Returning NA; talasso() will use mean(y) as forecast.")
        return(NA)
    } else if (p_selected > 1) {
        # talasso() always builds its second-step weights from the OLS fit on the
        # first-step-selected set (low-dimensional by construction). Tune lambda
        # against those same OLS-based weights so the selected penalty matches the
        # estimator that is actually deployed.
        best_tune <- train_lasso(xx, y, ada = TRUE, gamma = gamma, model = "ols",
                                 intercept = intercept, scale_x = scale_x,
                                 lambda_seq = lambda_seq, train_method = train_method,
                                 nlambda = nlambda, lambda_min_ratio = lambda_min_ratio,
                                 k = k, initial_window = initial_window,
                                 horizon = horizon, fixed_window = fixed_window,
                                 skip = skip)
        return(best_tune)
    } else {
        # Single-survivor path: weights also come from the OLS fit, matching the
        # second step of talasso() (init_est() with one column always picks OLS).
        init_result <- init_est(xx, y, gamma = gamma, intercept = intercept, scale_x = scale_x, model = "ols")
        w <- init_result$w
        coef_init <- init_result$coef_init

        if (is.null(lambda_seq)) {
            # lambda_max must match the corrected single-predictor soft-threshold
            # in lasso_weight_single() (which divides by the Gram term v). The
            # entry threshold that zeroes the coefficient is then the centered
            # cross-moment on the operating scale times max|b_init_std|^gamma,
            # mirroring the general adaptive grid in train_lasso().
            xv <- as.numeric(xx)
            xop <- if (scale_x) xv / sd_n(xv) else xv
            cross <- if (intercept) {
                abs(sum((xop - mean(xop)) * (y - mean(y)))) / n
            } else {
                abs(sum(xop * y)) / n
            }
            coef_slope <- if (intercept) coef_init[-1] else coef_init
            if (scale_x) coef_slope <- coef_slope * sd_n(xv)
            coef_max <- max(abs(coef_slope)^gamma)
            lambda_max <- coef_max * cross
            if (lambda_max == 0) lambda_max <- cross
            lambda_seq <- get_lambda_seq(lambda_max, lambda_min_ratio = lambda_min_ratio, nlambda = nlambda)
        }

        if (train_method %in% c("cv", "timeslice")) {
            if (train_method == "cv") {
                data_split <- list(train = list(), test = list())
                ind_seq <- 1:n
                seq_interval <- split(ind_seq, ceiling(seq_along(1:n) / (n / k)))
                for (j in 1:k) {
                    data_split$train[[j]] <- ind_seq[-seq_interval[[j]]]
                    data_split$test[[j]] <- ind_seq[seq_interval[[j]]]
                }
            } else if (train_method == "timeslice") {
                data_split <- create_time_slices(n, initial_window, horizon,
                                                    fixed_window, skip)
            }

            MSE <- rep(0, length(lambda_seq))
            for (i in seq_along(lambda_seq)) {
                lambda <- lambda_seq[i]
                for (j in seq_along(data_split$train)) {
                    y_j <- y[data_split$train[[j]]]
                    x_j <- xx[data_split$train[[j]], ]
                    y_p <- as.matrix(y[data_split$test[[j]]])
                    x_p <- xx[data_split$test[[j]], ]

                    result <- lasso_weight_single(x_j, y_j, lambda = lambda,
                                                  w = w, intercept = intercept, scale_x = scale_x)
                    a_ada <- as.numeric(result$ahat)
                    b_ada <- as.numeric(result$bhat)

                    if (intercept) {
                        coef_ada <- c(a_ada, b_ada)
                        mse_j <- colMeans((y_p - cbind(1, x_p) %*% coef_ada)^2)
                    } else {
                        coef_ada <- b_ada
                        mse_j <- mean((y_p - as.matrix(x_p) * coef_ada)^2)
                    }
                    MSE[i] <- MSE[i] + mse_j
                }
            }

            return(lambda_seq[which.min(MSE)])
        } else if (train_method %in% c("aic", "bic", "aicc", "hqc")) {
            Coef_hat <- matrix(0, p_selected + intercept, length(lambda_seq))
            Df <- rep(0, length(lambda_seq))

            for (ll in seq_along(lambda_seq)) {
                lambda <- lambda_seq[ll]
                result <- lasso_weight_single(xx, y, lambda = lambda,
                                              w = w, intercept = intercept, scale_x = scale_x)
                a_ada <- as.numeric(result$ahat)
                b_ada <- as.numeric(result$bhat)
                coef_ada <- if (intercept) c(a_ada, b_ada) else b_ada
                Coef_hat[, ll] <- coef_ada
                Df[ll] <- sum(b_ada != 0)
            }

            y_hat <- cbind(1, xx) %*% Coef_hat
            e_hat <- matrix(y, n, length(lambda_seq)) - y_hat
            mse <- colMeans(e_hat^2)

            nvar <- Df + intercept
            ic <- switch(train_method,
                "aic" = n * log(mse) + 2 * nvar,
                "bic" = n * log(mse) + log(n) * nvar,
                "aicc" = n * log(mse) + 2 * nvar + 2 * nvar * (nvar + 1) / (n - nvar - 1),
                "hqc" = n * log(mse) + 2 * nvar * log(log(n))
            )

            # See .train_ic(): drop saturated lambdas (mse = 0 or AICc divisor
            # <= 0) that would otherwise win which.min via a -Inf/NaN criterion.
            ic[!is.finite(ic)] <- Inf
            return(lambda_seq[which.min(ic)])
        } else {
            stop("Invalid train_method input.")
        }
    }
}
