#' Implement Adaptive Lasso via glmnet
#'
#' Estimation function. Tuning parameter inputs needed.\cr
#' Incorporates both high-dim(Lasso as initial estimator) and low-dim (OLS as initial estimator).
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param lambda Shrinkage tuning parameter for the adaptive step
#' @param lambda_init Shrinkage tuning parameter for the initial step Lasso estimation (If applicable)
#' @param gamma Parameter controlling the inverse of first step estimate
#' @param intercept A boolean: include an intercept term or not
#' @param scale_x A boolean: standardize the design matrix or not
#' @param model A string: model in the initial estimation, "lasso", "ridge", or "ols". Default is "lasso" for high-dimensional case and "ols" for low-dimensional case.
#'
#' @return A list contains estimated intercept and slope
#' \item{ahat}{Estimated intercept}
#' \item{bhat}{Estimated slope}
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' x <- matrix(rnorm(50 * 5), 50, 5)
#' y <- rnorm(50)
#' adalasso(x, y, lambda = 0.1)
#' }
adalasso <- function(
  x, y, lambda, lambda_init = NULL, gamma = 1, intercept = TRUE,
  scale_x = FALSE, model = NULL
) {
    n <- nrow(x)
    p <- ncol(x)

    init <- init_est(x, y, lambda_init, gamma, intercept, scale_x, model)
    w <- init$w
    model <- init$model

    # Adaptive Lasso estimation
    result <- glmnet::glmnet(x, y,
        lambda = lambda * sum(w) / p, penalty.factor = w,
        intercept = intercept, standardize = scale_x
    )
    ahat <- as.numeric(result$a0) # If intercept = FALSE, this will be 0
    bhat <- as.numeric(result$beta)

    return(list(ahat = ahat, bhat = bhat, model = model))
}

init_est <- function(x, y, lambda_init = NULL, gamma = 1, intercept = TRUE, scale_x = FALSE, model = NULL) {
    n <- nrow(x)
    p <- ncol(x)

    # Decide which model to use
    if (is.null(model)) {
        model <- if (p > 2 * sqrt(n)) "lasso" else "ols"
    }

    # Validate model choice
    if (!model %in% c("ols", "lasso", "ridge")) {
        stop("Error: model must be one of NULL, 'ols', 'lasso', or 'ridge'.")
    }

    # OLS not feasible in high dimension
    if (model == "ols" && p > n) {
        stop("Error: p > n, OLS is not feasible. Use Lasso or ridge instead.")
    }

    # Tune lambda_init for Lasso or ridge if not provided
    if (model %in% c("lasso", "ridge") && is.null(lambda_init)) {
        alpha <- ifelse(model == "lasso", 1, 0)
        lambda_init <- train_lasso(x, y, ada = FALSE, intercept = intercept, scale_x = scale_x, alpha = alpha, train_method = "cv")
    }

    # Fit chosen model
    if (model == "ols") {
        coef_init <- lsfit(x, y, intercept = intercept)$coef

        if (intercept) {
            b_temp <- coef_init[-1]
        } else {
            b_temp <- coef_init
        }
        b_temp[is.na(b_temp)] <- 1e-8
    } else if (model == "lasso") {
        lasso_result <- glmnet::glmnet(
            x, y,
            lambda = lambda_init,
            intercept = intercept,
            standardize = scale_x
        )

        b_temp <- as.numeric(lasso_result$beta)

        if (intercept) {
            coef_init <- c(as.numeric(lasso_result$a0), b_temp)
        } else {
            coef_init <- b_temp
        }
    } else if (model == "ridge") {
        ridge_result <- glmnet::glmnet(
            x, y,
            lambda = lambda_init,
            alpha = 0,
            intercept = intercept,
            standardize = scale_x
        )

        b_temp <- as.numeric(ridge_result$beta)

        if (intercept) {
            coef_init <- c(as.numeric(ridge_result$a0), b_temp)
        } else {
            coef_init <- b_temp
        }
    }

    # Adaptive weights live on the penalty scale. When scale_x = TRUE the
    # second step penalizes standardized coefficients, so the initial estimate
    # must be read on the standardized scale (theta = b * sd) for the weighting
    # to reflect importance rather than the predictors' units (scale invariance).
    b_std <- if (scale_x) b_temp * apply(x, 2, sd_n) else b_temp
    # Floor exact zeros (Lasso-screened variables) on the penalty scale, not the
    # original scale, so that their effectively-infinite weight is unit-free.
    b_std[b_std == 0] <- 1e-08
    w <- 1 / (abs(b_std)^gamma)

    return(list(w = w, coef_init = coef_init, model = model))
}
