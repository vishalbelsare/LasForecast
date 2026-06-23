#' Implement repeated Lasso estimation
#'
#' Estimation function. Tuning parameter inputs needed. First step adalasso estimate needed.\cr
#'
#' @param b_first First step adaptive lasso estimates
#' @param gamma Parameter controlling the inverse of first step estimate
#' @inheritParams lasso_weight_single
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
#' b_first <- rnorm(5)
#' talasso(x, y, b_first, lambda = 0.1)
#' }
talasso <- function(x,
                    y,
                    b_first,
                    lambda,
                    gamma = 1,
                    intercept = TRUE,
                    scale_x = FALSE) {
    p_full <- ncol(x)

    if (sum(b_first == 0) == p_full) {
        return(list(ahat = mean(y), bhat = rep(0, p_full)))
    }

    b_temp <- b_first

    xx <- x[, b_first != 0]
    coef_ols_second <- lsfit(xx, y, intercept = intercept)$coef

    if (intercept) {
        b_ols_second <- coef_ols_second[-1]
    } else {
        b_ols_second <- coef_ols_second
    }

    # As in init_est(): when scale_x = TRUE the penalty below acts on
    # standardized coefficients, so the OLS estimate must be read on the
    # standardized scale (theta = b * sd) for the weights to be unit-free.
    b_std <- if (scale_x) b_ols_second * apply(as.matrix(xx), 2, sd_n) else b_ols_second
    w <- 1 / (abs(b_std)^gamma)

    # Second Adaptive Lasso estimation

    if (sum(b_first != 0) == 1) {
        result <- lasso_weight_single(xx, y,
            lambda = lambda,
            w = w, intercept = intercept, scale_x = scale_x
        )
        ahat_second <- as.numeric(result$ahat)
        bhat_second <- as.numeric(result$bhat)
    } else {
        p_sel <- ncol(xx)
        result <- glmnet::glmnet(as.matrix(xx), y,
            lambda = lambda * sum(w) / p_sel,
            penalty.factor = w, intercept = intercept, standardize = scale_x
        )
        ahat_second <- as.numeric(result$a0)
        bhat_second <- as.numeric(result$beta)
    }

    b_temp[b_first != 0] <- bhat_second
    return(list(ahat = ahat_second, bhat = b_temp))
}
