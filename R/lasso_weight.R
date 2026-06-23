#' Implement weighted Lasso estimation (second step of adaptive lasso) via coordinate descent
#' for the case in which only one predictor remains
#'
#' Estimation function. Tuning parameter inputs needed.\cr
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param lambda Shrinkage tuning parameter
#' @param w weights
#' @param intercept A boolean: include an intercept term or not
#' @param scale_x A boolean: standardize the design matrix or not
#'
#' @return A list contains estimated intercept and slope
#' \item{ahat}{Estimated intercept}
#' \item{bhat}{Estimated slope}
lasso_weight_single <- function(x,
                                y,
                                lambda,
                                w = NULL,
                                intercept = TRUE,
                                scale_x = FALSE) {

    if(scale_x) {
        s_x <- sd_n(x)
        x <- x / s_x
    }
    if(is.null(w)) w = 1

    # Gram term for the single predictor on the operating scale. The weighted
    # soft-threshold acts on the OLS slope and must be divided by this (glmnet
    # divides each coordinate update by x_j'x_j / n); omitting it would solve
    # the problem for penalty lambda * w * v instead of lambda * w. With an
    # intercept the relevant moment is centered; without one it is uncentered.
    # When scale_x = TRUE the column was rescaled to unit centered variance, so
    # v = 1 (and the threshold is unchanged).
    v <- if (intercept) sd_n(x)^2 else mean(x^2)

    if(intercept) b_ols <- lsfit(x,y)$coefficients[-1]
    else b_ols <- lsfit(x,y,intercept = intercept)$coefficients

    bhat <- sign(b_ols) * pos(abs(b_ols) - lambda * w / v)
    if (intercept) ahat <- mean(y - x*bhat)
    else ahat = NULL

    # bhat above is on the standardized scale; map it back to the original
    # scale so callers can apply it to unscaled predictors.
    if(scale_x) bhat <- bhat / s_x

    return(list(ahat = ahat, bhat = bhat))

}

