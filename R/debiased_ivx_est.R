# Main function for debiased IVX estimator

#' Debiased IVX for predictive regression
#'
#' \deqn{y_t = w_{t-1} theta + u_t} where data is already aligned to incorporate the lagged regressor
#'
#' @param w Matrix of all regressors
#' @param y Vector of dependent variable
#' @param d_ind Index for inference targets
#' @param intercept Whether to include intercept in Lasso regression
#' @param standardize Whether to standardize the variables
#' @param c_z Parameter in constructing IV (Phillips and Lee, 2016) \deqn{z = \sum_{j=0}^{n-1} (1 - c_z / n^a)^j \Delta d_{t-j}}
#' @param a Parameter in constructing IV (Phillips and Lee, 2016)
#' @param standardize_iv Whether to standardize the IV
#' @param lambda_choice Choice of lambda for Lasso regression; List of length length(d_ind) + 1: each element = NULL or = a number if user has a specific choice of tuning parameter
#' @param lambda_seq pre-specified sequence of tuning parameter for parameter tuning; Useful in calibration of tuning parameter based on the rate conditions in the asymptotic theory; List of length length(d_ind) + 1: Each element = NULL or a vector of tuning parameters
#' @param train_method The parameter tuning method
#'     \itemize{
#'      \item \code{"timeslice"}: Time-slice cross-validation.
#'          By combining initial window, horizon, fixed window and skip, we can control the sample splitting.
#'          Roll_block: Setting initial_window = horizon = floor(nrow(x) / k), fixed_window = False, and skip = floor(nrow(x) / k) - 1.
#'          Period-by-period rolling: skip = 0.
#'      \item \code{"cv"}: Cross-validation based on block splits.
#'      \item \code{"cv_random"}: Cross-validation based on random splits.
#'      \item \code{"aic"}, \code{"bic"}, \code{"aicc"}, \code{"hqc"}: based on information criterion.
#'      }
#' @param nlambda number of candidate lambdas
#' @param lambda_min_ratio # lambda_min_ratio * lambda_max = lambda_min (default: 0.0001): Determines the search range of lambda
#' @param k k-fold cv if "cv" is chosen (default: 10)
#' @param initial_window length of initial window for "timeslice" method
#' @param horizon length of horizon for "timeslice" method
#' @param fixed_window whether to use fixed window for "timeslice" method
#' @param skip length of skip for "timeslice" method
#' @param joint_test boolean indicating whether to conduct a joint Wald test
#' @param R_mat restriction matrix for joint test
#' @param q_vec restriction vector for joint test
#' @param se_type type of standard error estimation
#'     \itemize{
#'      \item "iid": iid standard error
#'      \item "robust": robust standard error
#'      \item "HAC": HAC standard error
#'      }
#' @return A list contains
#' \item{theta_hat_las}{Estimate of delta from Lasso regression}
#' \item{theta_hat_ivx}{Estimate of delta from debiased IVX}
#' \item{sigma_hat}{Estimate of standard error of theta_hat_ivx}
#' \item{lambda_hat}{Chosen tuning parameter for Lasso regression}
#' \item{phi_hat}{Estimate of the frequency of 0s and L1 norm of std/nonstd coefficients in the second stage}
#'
#' @export
#'


debias_ivx <- function(
  w,
  y,
  d_ind,
  intercept = FALSE,
  standardize = TRUE,
  c_z = 5,
  a = 0.5,
  standardize_iv = TRUE,
  lambda_choice = vector("list", length(d_ind) + 1),
  lambda_seq = vector("list", length(d_ind) + 1),
  train_method = "timeslice",
  nlambda = 100,
  lambda_min_ratio = 0.0001,
  k = 10,
  initial_window = ceiling(nrow(w) * 0.7),
  horizon = 1,
  fixed_window = TRUE,
  skip = 0,
  se_type = "iid",
  joint_test = FALSE,
  R_mat = diag(length(d_ind)),
  q_vec = rep(0, length(d_ind))
) {
    n <- length(y)
    p_focal <- length(d_ind)
    p <- ncol(w)

    fit_lasso_args <- list(
        intercept = intercept,
        standardize = standardize,
        train_method = train_method,
        nlambda = nlambda,
        lambda_min_ratio = lambda_min_ratio,
        k = k,
        initial_window = initial_window,
        horizon = horizon,
        fixed_window = fixed_window,
        skip = skip
    )

    # Container for chosen tuning parameters
    lambda_hat <- rep(NA, length(d_ind) + 1)

    # ---- Step 1: Lasso Regression y on w ------
    lasso_result <- do.call(
        fit_lasso,
        c(list(w = w, y = y, lambda_choice = lambda_choice[[1]], lambda_seq = lambda_seq[[1]]), fit_lasso_args)
    )
    b_hat_las <- lasso_result$beta
    u_hat <- as.numeric(lasso_result$u)
    theta_hat_las <- b_hat_las[d_ind]
    lambda_hat[1] <- lasso_result$lambda

    # Make this function also a wrapper for Lasso regression
    if (is.null(d_ind)) {
        return(lasso_result)
    }

    # --------------------------------------------

    # ---- Step 2: IVX ----------------------------
    theta_hat_ivx <- rep(NA, p_focal)
    sigma_hat_ivx <- rep(NA, p_focal)

    # Container for second stage estimated coefficients
    # Three rows: frequency of 0s, L1 norm of std/nonstd.
    phi_hat <- matrix(NA, 3, p_focal)

    w_joint <- w[-1, d_ind, drop = FALSE]
    r_joint <- matrix(NA, nrow(w_joint), p_focal)

    for (i in 1:p_focal) {
        d <- w[, d_ind[i]]


        z <- generate_iv(d, n, a = a, c_z = c_z)
        # Normalize the IV
        if (standardize_iv) {
            z <- z / sd_n(z)
        }

        w_z <- w[-1, -d_ind[i], drop = FALSE]
        lasso_result <- do.call(
            fit_lasso,
            c(list(w = w_z, y = z, lambda_choice = lambda_choice[[i + 1]], lambda_seq = lambda_seq[[i + 1]]), fit_lasso_args)
        )
        b_hat_las_z <- as.numeric(lasso_result$beta)
        r_hat <- as.numeric(lasso_result$u)
        r_joint[, i] <- r_hat
        lambda_hat[i + 1] <- as.numeric(lasso_result$lambda)

        # glmnet reports the coefficients beta_j instead of beta_j * sd_j
        phi_hat[, i] <- c(
            mean(b_hat_las_z == 0),
            sum(abs(b_hat_las_z)),
            sum(abs(b_hat_las_z * apply(w_z, 2, sd_n)))
        )

        # # Generate debiased estimates
        if (se_type == "iid") {
            theta_hat_ivx[i] <- theta_hat_las[i] + (sum(r_hat * u_hat[-1])) / sum(r_hat * d[-1])
            omega_uu <- mean(u_hat^2)
            sigma_hat_ivx[i] <- sqrt(
                (omega_uu * sum(r_hat^2)) / (sum(r_hat * d[-1])^2)
            )
        } else if (se_type == "robust") {
            theta_hat_ivx[i] <- theta_hat_las[i] + (sum(r_hat * u_hat[-1])) / sum(r_hat * d[-1])
            sigma_hat_ivx[i] <- sqrt(
                sum((r_hat * u_hat[-1])^2) / (sum(r_hat * d[-1])^2)
            )
        } else {
            lrcov_du <- lrcov_est(u_hat[-1], diff(d), type = 1) # one-sided long-run covariance
            theta_hat_ivx[i] <- theta_hat_las[i] + (sum(r_hat * u_hat[-1]) - (n * lrcov_du)) / sum(r_hat * d[-1])
            omega_uu <- lrcov_est(u_hat, type = 0) # long-run covariance
            sigma_hat_ivx[i] <- sqrt(
                (omega_uu * sum(r_hat^2)) / (sum(r_hat * d[-1])^2)
            )
        }
    }

    if (joint_test == TRUE) {
        # The joint covariance must follow the same convention as the marginal
        # SE, so that for a single restriction the Wald statistic equals the
        # squared t-statistic. Under se_type = "iid" (default) this is the
        # homoskedastic form of eq. (2.19) in Gao et al. (2026),
        # Omega_{j,k} = sigma_u^2 sum(r_j r_k) / (sum(r_j w_j) sum(r_k w_k));
        # "robust" uses the heteroskedasticity-consistent sandwich; the HAC
        # branch uses the long-run variance of u (matching its marginal SE).
        omega_uu_joint <- if (se_type == "iid") {
            mean(u_hat^2)
        } else if (se_type == "robust") {
            NULL
        } else {
            lrcov_est(u_hat, type = 0)
        }
        cov_matrix <- matrix(NA, p_focal, p_focal)
        for (i in 1:p_focal) {
            for (j in 1:p_focal) {
                num <- if (se_type == "robust") {
                    sum(r_joint[, i] * r_joint[, j] * (u_hat[-1]^2))
                } else {
                    omega_uu_joint * sum(r_joint[, i] * r_joint[, j])
                }
                cov_matrix[i, j] <- num /
                    (sum(r_joint[, i] * w_joint[, i]) * sum(r_joint[, j] * w_joint[, j]))
            }
        }
        test_vec <- R_mat %*% as.matrix(theta_hat_ivx, ncol = 1) - q_vec
        R_cov <- R_mat %*% cov_matrix %*% t(R_mat)
        wald_stat <- as.numeric(t(test_vec) %*% solve(R_cov) %*% test_vec)
        p_value_wald <- pchisq(wald_stat, df = nrow(R_mat), lower.tail = FALSE)
    } else {
        wald_stat <- NA
        p_value_wald <- NA
    }

    # --------------------------------------------
    output_list <- list(
        theta_hat_las = theta_hat_las,
        theta_hat_ivx = theta_hat_ivx,
        sigma_hat_ivx = sigma_hat_ivx,
        lambda_hat = lambda_hat,
        phi_hat = phi_hat,
        b_hat_las = b_hat_las,
        u_hat = u_hat,
        wald_stat = wald_stat,
        p_value_wald = p_value_wald
    )

    return(output_list)
}


#' Self-generated IVs
#'
#' @param d Vector of regressor values
#' @param n Sample size
#' @param a Parameter in constructing IV
#' @param c_z Parameter in constructing IV
#'
#' @return Instrumental variable vector of length n-1
#'
#' @keywords internal
generate_iv <- function(d, n, a, c_z = 5) {
    delta_d <- c(0, diff(d))
    d_mat <- toeplitz(delta_d)
    d_mat[upper.tri(d_mat)] <- 0
    const_mat <- (1 - (c_z / n^a))^matrix(0:(n - 1), n, n, byrow = TRUE)
    const_mat[upper.tri(const_mat)] <- 0
    z <- rowSums(const_mat * d_mat)[-1] # Dimension of Z is n - 1

    return(z)
}


#' Run Lasso estimation
#'
#' @param w Matrix of all regressors
#' @param y Vector of dependent variable
#' @param intercept Whether to include intercept in Lasso regression
#' @param standardize Whether to standardize the variables
#' @inheritParams debias_ivx
#'
#' @return coefficients of Lasso regression and residuals
#'
#' @keywords internal
#'

# Rewrite the function to avoid repeition of function calling.
fit_lasso <- function(
  w, y,
  intercept = FALSE,
  standardize = TRUE,
  lambda_choice = NULL,
  lambda_seq = NULL,
  train_method = "timeslice",
  nlambda = 100,
  lambda_min_ratio = 0.0001,
  k = 10,
  initial_window = ceiling(nrow(w) * 0.7),
  horizon = 1,
  fixed_window = TRUE,
  skip = 0
) {
    if (is.null(lambda_choice)) {
        train_arg <- list(
            x = w,
            y = y,
            ada = FALSE,
            intercept = intercept,
            scale_x = standardize,
            lambda_seq = lambda_seq,
            train_method = train_method,
            nlambda = nlambda,
            lambda_min_ratio = lambda_min_ratio,
            k = k,
            initial_window = initial_window,
            horizon = horizon,
            fixed_window = fixed_window,
            skip = skip
        )
        lambda_lasso <- do.call(train_lasso, train_arg)
    } else {
        lambda_lasso <- lambda_choice
    }
    result <- glmnet::glmnet(w,
        y,
        lambda = lambda_lasso,
        intercept = intercept,
        standardize = standardize
    )
    b_hat_las <- result$beta
    u_hat <- y - w %*% b_hat_las - result$a0 # result$a0 = 0 if intercept = FALSE

    return(
        list(beta = b_hat_las, u = u_hat, lambda = lambda_lasso)
    )
}

