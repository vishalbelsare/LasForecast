# ---- Generate lambda Sequence as in glmnet ----
# Details refer to Friedman, Hastie and Tibshirani (2010, Sec. 2.5.) Journal of Statistical Software
# and the online posts:
# https://stackoverflow.com/questions/23686067/default-lambda-sequence-in-glmnet-for-cross-validation
# https://stackoverflow.com/questions/25257780/how-does-glmnet-compute-the-maximal-lambda-value/


#' standard deviation of x
#' @param x A numeric vector
#' @return SD of x
sd_n <- function(x) {
    sqrt(sum((x - mean(x))^2) / length(x))
}


#' Get the maximum lambda value for lasso
#'
#' @param x Predictor matrix
#' @param y Response vector
#' @param scale_x Boolean to scale x
#' @param intercept Boolean; if TRUE (and \code{scale_x = FALSE}) center x and y
#'   so the grid top matches glmnet's centered lambda_max.
#' @param nlambda Number of lambdas
#' @param lambda_min_ratio Ratio of min/max lambda
#' @param per_coord Boolean; if TRUE return the per-coordinate entry thresholds
#'   \eqn{|x_j'(y - \bar y)| / n} (a length-\code{ncol(x)} vector) instead of
#'   their maximum. The adaptive grid in \code{train_lasso} weights each
#'   coordinate before taking the max, so it needs the vector.
#'
#' @return Scalar \code{lambda_max} (the maximum), or the per-coordinate vector
#'   when \code{per_coord = TRUE}.
#'
#' @keywords internal
#'

get_lasso_lambda_max <- function(x, y, scale_x = FALSE, intercept = TRUE, nlambda = 100, lambda_min_ratio = 0.0001, per_coord = FALSE) {
    n <- nrow(x)
    p <- ncol(x)
    y <- as.numeric(y)

    if (scale_x) {
        sx <- scale(x, scale = apply(x, 2, sd_n))
        sx <- as.matrix(sx, ncol = p, nrow = n)

        cross <- abs(colSums(sx * y)) / n
    } else {
        # With an intercept glmnet evaluates the entry gradient at the centered
        # design (a = ybar, b = 0); center x and y so the grid top matches its
        # lambda_max instead of the inflated uncentered cross-product.
        if (intercept) {
            x <- scale(x, center = TRUE, scale = FALSE)
            y <- y - mean(y)
        }
        cross <- abs(colSums(x * y)) / n
    }

    # cross_j is coordinate j's entry threshold at the all-zero solution; the
    # plain-Lasso grid top is their max. per_coord exposes the vector so the
    # adaptive grid can apply the weights 1 / |theta_init_j|^gamma first.
    if (per_coord) cross else max(cross)
}


#' Get the lambda sequence
#'
#' @param lambda_max Max lambda
#' @param lambda_min_ratio Ratio
#' @param nlambda Number of lambdas
#'
#' @return lambda_seq
#'
#' @keywords internal
#'

get_lambda_seq <- function(lambda_max, lambda_min_ratio = 0.0001, nlambda = 100) {
    lambda_min <- lambda_min_ratio * lambda_max
    lambda_seq <- exp(seq(log(lambda_max), log(lambda_min), length.out = nlambda))

    return(lambda_seq)
}

# Positive-part operator (x)_+ = max(x, 0); vectorized for safety.
pos <- function(x) pmax(x, 0)
# ---- end of lambda sequence ----
