#' post-selection estimation
#'
#' @param x Predictor matrix (n-by-p matrix)
#' @param y Response variable
#' @param coef_est shrinkage estimation results
#' @param intercept A boolean: include an intercept term or not
#' @param scale_x A boolean: standardize the design matrix or not
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
#' coef_est <- c(0, 1, 0, 0, -0.5, 0)
#' post_lasso(x, y, coef_est)
#' }
#'
post_lasso <- function(x, y, coef_est, intercept = TRUE, scale_x = FALSE) {

    p <- ncol(x)
    if(scale_x) x = scale(x, center = FALSE, scale = apply(x, 2, sd_n) )

    if (intercept)
        b_est <- coef_est[-1]
    else
        b_est <- coef_est

    if (sum(b_est != 0) == 0)
        return(coef_est)
    else
        coef_post_temp <- lsfit(x[, b_est != 0], y, intercept = intercept)$coef

    if (intercept) {
        b_post_temp <- coef_post_temp[-1]
        b_post <- rep(0, p)
        b_post[b_est != 0] <- b_post_temp
        coef_post <- c(coef_post_temp[1], b_post)
    } else {
        coef_post <- rep(0, p)
        coef_post[b_est != 0] <- coef_post_temp
    }

    return(coef_post)

}
