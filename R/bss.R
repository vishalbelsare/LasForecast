# Best subset selection via mixed integer optimization (MIO), following
# Bertsimas, King, and Mazumder (2016, Annals of Statistics).
#
# Requires the 'gurobi' R package, which is not on CRAN: it ships with the
# Gurobi Optimizer (https://www.gurobi.com) and is installed from the local
# Gurobi distribution. Declared in Suggests and checked at run time.

#' Best Subset Selection for Predictive Regression
#'
#' Computes the best subset selection estimator
#' \deqn{\min_\beta \frac{1}{2}\|y - X\beta\|_2^2 \quad \text{s.t.} \quad
#'   \|\beta\|_0 \le k}
#' by mixed integer optimization (MIO) with the Gurobi solver, following
#' Bertsimas, King, and Mazumder (2016). For each subset size the MIO is
#' warm-started by a discrete first-order algorithm (their Algorithm 2 with
#' polishing). The overdetermined case \eqn{p < n} uses their formulation
#' (2.5) (quadratic in \eqn{p} variables); the high-dimensional case
#' \eqn{p \ge n} uses formulation (2.6) (quadratic in \eqn{n} variables).
#'
#' When \code{k} is \code{NULL} (default), the subset size is selected by
#' information criterion over \eqn{k = 1, \dots, k_{max}} (and \eqn{k = 0},
#' the intercept-only model, when \code{include_k0 = TRUE}). The intercept is
#' handled exactly by centering (Frisch--Waugh): the MIO runs on demeaned
#' data and the intercept is recovered as
#' \eqn{\hat\alpha = \bar y - \bar x'\hat\beta}; it is not counted in
#' \eqn{k}.
#'
#' The problem is solved on the raw data throughout (centered when an
#' intercept is included): no standardization or rescaling of the
#' predictors is applied at any stage, so the estimator inherits the exact
#' scale-equivariance of least squares on the selected subset. The MIO
#' imposes no heuristic coefficient bounds -- the SOS-1 formulation needs
#' no big-M constants (their formulation (2.4)) -- so the solver returns
#' the exact constrained least squares optimum (up to solver tolerances)
#' in every regime. The formulation is strengthened only by provable
#' constraints derived from the warm-start objective value \eqn{UB}:
#' per-coefficient bounds from the least squares level set when
#' \eqn{p < n} and \eqn{X'X} is invertible (Section 2.3.2 of the paper,
#' closed form for least squares), and fitted-value bounds
#' \eqn{|\eta_i - y_i| \le \sqrt{2 UB}} when \eqn{p \ge n}. Predictors
#' that are constant over the estimation window cannot enter the fit and
#' receive a zero coefficient.
#'
#' @param x Predictor data: an \eqn{n \times p}{n x p} numeric matrix,
#'   data.frame, tibble, or tsibble. For data.frame/tibble/tsibble inputs,
#'   \code{y} can be a column name and date columns are auto-detected.
#' @param y Response vector (length \eqn{n}), or a column name string
#'   when \code{x} is a data.frame/tibble/tsibble.
#' @param k Candidate subset size(s). One of: \code{NULL} (default) --
#'   \eqn{k} is selected by \code{ic_type} over \eqn{0/1, \dots, k_{max}};
#'   an integer scalar in \code{0:p} -- the model with at most \eqn{k}
#'   active predictors is estimated directly (\code{k = 0} is the null
#'   model: the prevailing mean when \code{intercept = TRUE}, the zero
#'   forecast otherwise); or an integer vector (a subset of \code{0:p}) --
#'   \eqn{k} is selected by \code{ic_type} within that candidate set.
#' @param ic_type Information criterion for subset size selection:
#'   \code{"bic"} (default), \code{"aic"}, \code{"aicc"}, or \code{"hqc"}.
#' @param intercept Logical; include an (unpenalized) intercept
#'   (default \code{TRUE}).
#' @param include_k0 Logical; also consider \eqn{k = 0} (no predictors:
#'   the prevailing mean when \code{intercept = TRUE}, the zero forecast
#'   otherwise) in the information criterion comparison (default \code{TRUE}).
#'   Ignored when \code{k} is supplied.
#' @param k_max Largest subset size searched when \code{k} is \code{NULL}.
#'   Default (and upper cap) \code{min(p, n - 1 - intercept)}, so the least
#'   squares fit on the selected subset keeps positive residual degrees of
#'   freedom and the information criteria remain well defined.
#' @param predictive Logical. If \code{TRUE}, automatically lag-align the data
#'   so that row \eqn{t} pairs \eqn{y_t} with \eqn{x_{t-1}}. If \code{FALSE},
#'   the data is used as-is. Default \code{NULL} auto-detects: \code{TRUE}
#'   when date information is present, \code{FALSE} otherwise.
#' @param date_col Optional date column name (for data.frame/tibble input).
#' @param a Multiplier on the theory-guided bound radii (default 1). The
#'   MIO bounds are \eqn{|\beta_j| \le |\hat\beta_{OLS,j}| +
#'   a\sqrt{2(UB - g_{OLS})(X'X)^{-1}_{jj}}} when \eqn{p < n}, and
#'   \eqn{|\eta_i - y_i| \le a\sqrt{2 UB}},
#'   \eqn{\|\eta\|_1 \le \|y\|_1 + a\sqrt{2 n UB}} when \eqn{p \ge n}.
#'   \code{a = 1} gives the exact level-set bounds (the solution is provably
#'   never cut off); \code{a > 1} loosens them (still exact, possibly slower
#'   certification); \code{a < 1} tightens them beyond their provable level,
#'   which can speed up the solve but voids the exactness guarantee. No
#'   effect when the model is solved without finite bounds (singular
#'   \eqn{X'X}).
#' @param polish Logical; re-fit the first-order solution by least squares on
#'   its active set before passing it to the MIO (default \code{TRUE}),
#'   as recommended in Section 3.2 of the paper.
#' @param tol Convergence tolerance of the first-order algorithm on the
#'   objective (default \code{1e-4}, the value used in the paper).
#' @param max_iter Maximum number of first-order iterations (default 1000).
#' @param time_limit Gurobi time limit in seconds \emph{per MIO solve}
#'   (default 60). If the limit binds, the incumbent (best feasible) solution
#'   is used and a warning reports the remaining MIP gap. Proving optimality
#'   can be slow in the high-dimensional regime (\eqn{p \ge n}), where the
#'   exact formulation has weak relaxation bounds; the incumbent is
#'   typically optimal long before the certificate arrives (see Figure 2 of
#'   the paper).
#' @param verbose Logical; print Gurobi solver output (default \code{FALSE}).
#' @param ... Additional arguments (ignored). Accepted for interface
#'   compatibility with the other estimators.
#'
#' @return An object of class \code{"lasforecast_model"} (a list) with:
#'   \describe{
#'     \item{\code{fit}}{Estimation details: selected \code{k}, solver
#'       \code{status}, \code{mipgap}, \code{objval}, and (when \code{k} is
#'       selected by IC) the search path \code{ic_path}, a data.frame with
#'       one row per candidate \eqn{k}.}
#'     \item{\code{tuning_param}}{The selected subset size \eqn{\hat k}.}
#'     \item{\code{coefficients}}{Numeric vector of estimated coefficients
#'       (intercept first when \code{intercept = TRUE}).}
#'     \item{\code{method}}{\code{"bss"}.}
#'   }
#'   S3 methods: \code{\link{predict.lasforecast_model}},
#'   \code{\link[=print.lasforecast_model]{print}},
#'   \code{\link[=coef.lasforecast_model]{coef}}.
#'
#' @references
#' Bertsimas, D., King, A., and Mazumder, R. (2016).
#' "Best Subset Selection via a Modern Optimization Lens."
#' \emph{The Annals of Statistics}, 44(2), 813--852.
#'
#' @seealso \code{\link{lasso}} for Lasso-family estimation,
#'   \code{\link{roll_predict}} for rolling-window forecasting.
#'
#' @examples
#' \dontrun{
#' # Requires the Gurobi solver and its R package.
#' x <- matrix(rnorm(100 * 6), 100, 6)
#' y <- x[, 1] - 0.5 * x[, 2] + rnorm(100)
#'
#' # Subset size selected by BIC
#' mod <- bss(x, y)
#' coef(mod)
#'
#' # Fixed subset size
#' mod2 <- bss(x, y, k = 2)
#' }
#'
#' @export
bss <- function(x, y = NULL, k = NULL,
                ic_type = c("bic", "aic", "aicc", "hqc"),
                intercept = TRUE, include_k0 = TRUE, k_max = NULL,
                predictive = NULL, date_col = NULL, a = 1,
                polish = TRUE, tol = 1e-4, max_iter = 1000,
                time_limit = 60, verbose = FALSE, ...) {
    if (!requireNamespace("gurobi", quietly = TRUE)) {
        stop("bss() requires the 'gurobi' package, which ships with the ",
             "Gurobi Optimizer (https://www.gurobi.com). See the Gurobi ",
             "documentation for installing its R interface.")
    }
    ic_type <- match.arg(ic_type)
    if (length(a) != 1 || !is.numeric(a) || a <= 0) {
        stop("'a' must be a positive scalar.")
    }

    inp <- prepare_input(x, y, date_col)
    x <- inp$x_mat
    y <- inp$y_vec

    if (is.null(predictive)) predictive <- !is.null(inp$dates)
    if (predictive) {
        n_raw <- nrow(x)
        x <- x[-n_raw, , drop = FALSE]
        y <- y[-1]
    }

    n <- nrow(x)
    p <- ncol(x)

    # Center out the intercept (Frisch-Waugh): the slopes of the best subset
    # problem with a free intercept equal those of the no-intercept problem on
    # demeaned data, and alpha_hat = mean(y) - mean(x)' beta_hat.
    if (intercept) {
        x_means <- colMeans(x)
        y_mean <- mean(y)
        xc <- sweep(x, 2, x_means)
        yc <- y - y_mean
    } else {
        xc <- x
        yc <- y
    }

    # No standardization or rescaling anywhere: the problem is solved on the
    # raw (centered) data, so the estimator inherits the exact
    # scale-equivariance of least squares on the selected subset. Constant
    # columns (zero after centering) cannot enter the fit: their coefficient
    # is pinned to zero so no unidentified value leaks into forecasts.
    zero_var <- unname(colSums(xc^2)) == 0

    fit_k <- function(kk) {
        # k = p with p < n imposes no cardinality restriction: plain OLS.
        if (kk == p && p < n) {
            beta <- unname(suppressWarnings(
                stats::lsfit(xc, yc, intercept = FALSE)
            )$coef)
            beta[!is.finite(beta)] <- 0
            f <- list(beta = beta, status = "OPTIMAL", mipgap = 0,
                      objval = 0.5 * sum((yc - xc %*% beta)^2))
        } else {
            b0 <- bss_first_order(yc, xc, kk, tol = tol,
                                  max_iter = max_iter, polish = polish)
            f <- bss_mio(yc, xc, kk, b0, a = a,
                         time_limit = time_limit, verbose = verbose)
        }
        f$beta[zero_var] <- 0
        f
    }

    # Candidate set: a user scalar (fixed fit), a user vector (IC selection
    # within the set), or NULL (IC selection over 0/1, ..., k_max).
    if (!is.null(k)) {
        if (length(k) == 0 || !is.numeric(k) || anyNA(k) ||
                any(k != floor(k)) || any(k < 0) || any(k > p)) {
            stop("'k' must be NULL, an integer scalar in 0:ncol(x), or a ",
                 "vector of such values.")
        }
        ks <- sort(unique(as.integer(k)))
        select <- length(ks) > 1L
    } else {
        select <- TRUE
    }

    k_cap <- min(p, n - 1L - as.integer(intercept))
    if (is.null(k)) {
        if (is.null(k_max)) {
            k_max <- k_cap
        } else if (k_max > k_cap) {
            warning("'k_max' truncated to ", k_cap, " so that the subset ",
                    "least squares fit keeps positive residual degrees of ",
                    "freedom (information criteria are undefined at a ",
                    "perfect fit).")
            k_max <- k_cap
        }
        ks <- seq_len(k_max)
        if (include_k0) ks <- c(0L, ks)
    } else if (select && any(ks > k_cap)) {
        warning("Candidate values k > ", k_cap, " dropped: the subset ",
                "least squares fit must keep positive residual degrees of ",
                "freedom for the information criteria to be well defined.")
        ks <- ks[ks <= k_cap]
        if (length(ks) == 0) {
            stop("No candidate 'k' at or below ", k_cap, ".")
        }
    }

    fits <- lapply(ks, function(kk) {
        if (kk == 0L) {
            # Null model: no solver call needed.
            list(beta = rep(0, p), status = "OPTIMAL", mipgap = 0,
                 objval = 0.5 * sum(yc^2))
        } else {
            fit_k(kk)
        }
    })

    if (select) {
        mse <- vapply(fits, function(f) mean((yc - xc %*% f$beta)^2),
                      numeric(1))
        df <- vapply(fits, function(f) sum(f$beta != 0), numeric(1)) +
            as.integer(intercept)
        status <- vapply(fits, function(f) f$status, character(1))
        ic <- ic_eval(mse, df, n, ic_type)
        ic_path <- data.frame(k = ks, df = df, mse = mse, ic = ic,
                              status = status)
        sel <- which.min(ic)
    } else {
        ic_path <- NULL
        sel <- 1L
    }
    k_hat <- ks[sel]

    fit_sel <- fits[[sel]]
    not_opt <- vapply(fits, function(f) !identical(f$status, "OPTIMAL"),
                      logical(1))
    if (any(not_opt)) {
        warning("Gurobi did not certify optimality for k = ",
                paste(ks[not_opt], collapse = ", "),
                " (status: ", paste(unique(vapply(fits[not_opt],
                    function(f) f$status, character(1))), collapse = ", "),
                "). The best feasible solution found is used; consider ",
                "increasing 'time_limit'.")
    }

    beta <- fit_sel$beta
    coefs <- if (intercept) c(y_mean - sum(x_means * beta), beta) else beta

    structure(
        list(
            fit = list(k = k_hat, status = fit_sel$status,
                       mipgap = fit_sel$mipgap, objval = fit_sel$objval,
                       ic_type = if (select) ic_type else NULL,
                       ic_path = ic_path),
            tuning_param = k_hat,
            coefficients = coefs,
            method = "bss"
        ),
        class = "lasforecast_model"
    )
}


#' Discrete first-order algorithm for best subset selection
#'
#' Algorithm 2 of Bertsimas, King, and Mazumder (2016): projected gradient
#' steps with the hard-thresholding operator \eqn{H_k} (their Proposition 3)
#' and an exact line search (the objective is quadratic in the step size).
#' Following the paper, the algorithm exits with the best \emph{k-sparse}
#' iterate \eqn{\eta_{m^*}}, not the line-search combination (which can have
#' up to \eqn{2k} nonzeros), and the convergence criterion is
#' \eqn{|g(\eta_{m+1}) - g(\eta_m)| \le} \code{tol}. The algorithm is run
#' from two deterministic starting supports (raw-magnitude and scale-free
#' coefficient rankings; the paper itself uses multi-start, Section 5.2.1)
#' and the better run is kept. If \code{polish = TRUE}, the active set is
#' re-fitted by least squares (their Section 3.2).
#'
#' The data is assumed to be already centered when an intercept is intended.
#'
#' @param y Response vector.
#' @param x Predictor matrix.
#' @param k Subset size.
#' @param tol,max_iter Convergence control.
#' @param polish Logical; least squares re-fit on the active set.
#'
#' @return A k-sparse coefficient vector (length \code{ncol(x)}).
#'
#' @noRd
bss_first_order <- function(y, x, k, tol = 1e-4, max_iter = 1000,
                            polish = TRUE) {
    n <- nrow(x)
    p <- ncol(x)

    # Lipschitz constant of the gradient: largest eigenvalue of X'X,
    # computed from the smaller Gram matrix (X'X and XX' share it).
    gram <- if (p <= n) crossprod(x) else tcrossprod(x)
    L <- max(eigen(gram, symmetric = TRUE, only.values = TRUE)$values)

    # Base coefficients for initialization: unrestricted OLS when p < n,
    # marginal regression otherwise. Collinear or constant columns yield
    # NA/NaN; treat them as zero.
    if (p < n) {
        b_full <- unname(suppressWarnings(
            stats::lsfit(x, y, intercept = FALSE)
        )$coef)
    } else {
        b_full <- as.vector(crossprod(x, y)) / colSums(x^2)
    }
    b_full[!is.finite(b_full)] <- 0

    # Two deterministic starting supports, in the spirit of the paper's
    # multi-start practice (Section 5.2.1): threshold by raw coefficient
    # magnitude, and by the scale-free ranking |b_j| * ||X_j|| (marginal
    # correlation strength), which is robust when predictor scales are
    # heterogeneous. The rankings only pick starting supports; the data,
    # the iterations, and the estimator all stay on the raw scale.
    threshold_to <- function(rank_stat) {
        b <- b_full
        b[-order(rank_stat, decreasing = TRUE)[seq_len(k)]] <- 0
        b
    }
    inits <- unique(list(
        threshold_to(abs(b_full)),
        threshold_to(abs(b_full) * sqrt(colSums(x^2)))
    ))

    run_from <- function(b) {
        # The initialization is itself k-sparse: protect it in the
        # best-iterate tracking so a poor first projected-gradient step
        # cannot displace it.
        eta_best <- b
        g_eta_best <- 0.5 * sum((y - x %*% b)^2)
        g_eta_prev <- Inf

        for (m in seq_len(max_iter)) {
            r <- as.vector(y - x %*% b)
            # eta_m = H_k(b - grad/L): keep the k largest entries in
            # absolute value
            eta <- b + as.vector(crossprod(x, r)) / L
            eta[-order(abs(eta), decreasing = TRUE)[seq_len(k)]] <- 0

            g_eta <- 0.5 * sum((y - x %*% eta)^2)
            if (g_eta < g_eta_best) {
                g_eta_best <- g_eta
                eta_best <- eta
            }
            if (abs(g_eta_prev - g_eta) <= tol) break
            g_eta_prev <- g_eta

            # Exact line search: g(b + lambda * d) is quadratic in lambda.
            d <- eta - b
            xd <- as.vector(x %*% d)
            denom <- sum(xd^2)
            if (denom <= 0) break
            b <- b + (sum(r * xd) / denom) * d
        }
        list(eta = eta_best, g = g_eta_best)
    }

    runs <- lapply(inits, run_from)
    beta <- runs[[which.min(vapply(runs, `[[`, numeric(1), "g"))]]$eta

    if (polish && any(beta != 0)) {
        act <- beta != 0
        coef_act <- suppressWarnings(
            stats::lsfit(x[, act, drop = FALSE], y, intercept = FALSE)
        )$coef
        coef_act[!is.finite(coef_act)] <- 0
        beta[act] <- coef_act
    }
    beta
}


#' Best subset selection MIO for a given subset size
#'
#' Solves the cardinality-constrained least squares problem (no intercept;
#' center the data beforehand; raw scale throughout) with Gurobi. The
#' binary z here is the complement of the paper's: z_i = 1 forces
#' beta_i = 0 through an SOS-1 constraint on (beta_i, z_i), and
#' sum(z) >= p - k caps the active set at k. Because the SOS-1 formulation
#' needs no big-M constants (formulation (2.4) of Bertsimas, King, and
#' Mazumder, 2016), no heuristic coefficient bound is ever imposed: the
#' returned solution is the exact constrained least squares optimum up to
#' solver tolerances.
#'
#' The formulation is strengthened only by constraints that provably hold
#' for every candidate optimum, derived from the warm-start objective
#' UB = g(b0) (every optimum satisfies g(beta) <= UB):
#' \itemize{
#'   \item p < n with invertible X'X (their formulation (2.5), quadratic
#'     in p): the level set {g <= UB} is the ellipsoid around OLS with
#'     closed-form support function (their Section 2.3.2), giving
#'     per-coefficient bounds
#'     |beta_j| <= |beta_ols_j| + sqrt(2 (UB - g_ols) (X'X)^{-1}_jj) and
#'     ||beta||_1 <= sum of the k largest such bounds, linearized through
#'     the split beta = beta+ - beta- (no SOS-1 needed on the split: with
#'     a "<=" budget the minimal decomposition is always feasible).
#'   \item p < n with singular X'X: no valid finite box exists in the
#'     degenerate directions; the model is solved with free beta
#'     (formulation (2.4)).
#'   \item p >= n (their formulation (2.6), quadratic in n): free beta
#'     with eta = X beta and the provable fitted-value constraints
#'     |eta_i - y_i| <= sqrt(2 UB) and ||eta||_1 <= ||y||_1 +
#'     sqrt(2 n UB), the latter linearized through eta = eta+ - eta-.
#' }
#'
#' The multiplier \code{a} rescales the level-set radius terms in all of
#' these bounds (a = 1: provable; a > 1: looser but still valid; a < 1:
#' tighter than provable, exactness no longer guaranteed).
#'
#' @param y Centered response vector.
#' @param x Centered predictor matrix (raw scale; no standardization).
#' @param k Subset size.
#' @param b0 k-sparse warm start from \code{bss_first_order}.
#' @param a Multiplier on the bound radii.
#' @param time_limit Gurobi time limit (seconds).
#' @param verbose Logical; print solver output.
#'
#' @return List with \code{beta}, solver \code{status}, \code{mipgap},
#'   and \code{objval} (equal to \eqn{\frac{1}{2}\|y - X\beta\|^2}).
#'
#' @noRd
bss_mio <- function(y, x, k, b0, a = 1, time_limit = 60, verbose = FALSE) {
    n <- nrow(x)
    p <- ncol(x)

    # b0 = 0 happens iff X'y = 0, in which case beta = 0 is the global
    # optimum of the cardinality-constrained problem for every k.
    if (max(abs(b0)) == 0) {
        return(list(beta = rep(0, p), status = "OPTIMAL", mipgap = 0,
                    objval = 0.5 * sum(y^2)))
    }

    z0 <- as.numeric(b0 == 0)
    xty <- as.vector(crossprod(x, y))
    ub_obj <- 0.5 * sum((y - x %*% b0)^2)
    high_dim <- p >= n

    model <- list(modelsense = "min", objcon = 0.5 * sum(y^2))

    if (!high_dim) {
        # Provable per-coefficient bounds when X'X is invertible: {g <= UB}
        # is an ellipsoid around OLS with closed-form support function.
        ch <- tryCatch(chol(crossprod(x)), error = function(e) NULL)
        if (!is.null(ch)) {
            g_inv <- chol2inv(ch)
            beta_ols <- as.vector(g_inv %*% xty)
            g_ols <- 0.5 * sum((y - x %*% beta_ols)^2)
            r2 <- max(2 * (ub_obj - g_ols), 0)
            m_u_vec <- abs(beta_ols) + a * sqrt(r2 * pmax(diag(g_inv), 0))
            m_l <- sum(sort(m_u_vec, decreasing = TRUE)[seq_len(k)])

            # Formulation (2.5): variables (beta, z, beta+, beta-),
            # quadratic in p.
            # Rows: [1] sum(z) >= p - k; [2..p+1] beta - beta+ + beta- = 0;
            # [p+2] sum(beta+ + beta-) <= M_l.
            nv <- 4 * p
            model$obj <- c(-xty, rep(0, 3 * p))
            model$Q <- Matrix::bdiag(
                crossprod(x) / 2,
                Matrix::Matrix(0, 3 * p, 3 * p, sparse = TRUE)
            )
            model$A <- Matrix::sparseMatrix(
                i = c(rep(1, p), rep(2:(p + 1), 3), rep(p + 2, 2 * p)),
                j = c((p + 1):(2 * p),
                      1:p, (2 * p + 1):(3 * p), (3 * p + 1):(4 * p),
                      (2 * p + 1):(4 * p)),
                x = c(rep(1, p), rep(1, p), rep(-1, p), rep(1, p),
                      rep(1, 2 * p)),
                dims = c(p + 2, nv)
            )
            model$rhs <- c(p - k, rep(0, p), m_l)
            model$sense <- c(">", rep("=", p), "<")
            model$lb <- c(-m_u_vec, rep(0, 3 * p))
            model$ub <- c(m_u_vec, rep(1, p), m_u_vec, m_u_vec)
            model$vtype <- c(rep("C", p), rep("B", p), rep("C", 2 * p))
            model$start <- c(b0, z0, pmax(b0, 0), pmax(-b0, 0))
        } else {
            # Singular X'X: free beta, SOS-1 only (formulation (2.4)).
            nv <- 2 * p
            model$obj <- c(-xty, rep(0, p))
            model$Q <- Matrix::bdiag(
                crossprod(x) / 2,
                Matrix::Matrix(0, p, p, sparse = TRUE)
            )
            model$A <- Matrix::sparseMatrix(
                i = rep(1, p), j = (p + 1):(2 * p), x = rep(1, p),
                dims = c(1, nv)
            )
            model$rhs <- p - k
            model$sense <- ">"
            model$lb <- c(rep(-Inf, p), rep(0, p))
            model$ub <- c(rep(Inf, p), rep(1, p))
            model$vtype <- c(rep("C", p), rep("B", p))
            model$start <- c(b0, z0)
        }
    } else {
        # Formulation (2.6): variables (beta, z, eta, eta+, eta-) with
        # eta = X beta; quadratic in n; free beta. Provable eta bounds:
        # 0.5 ||y - eta||^2 <= UB gives the coordinate-wise box around y
        # and the l1 cap ||eta||_1 <= ||y||_1 + sqrt(2 n UB).
        nv <- 2 * p + 3 * n
        eta0 <- as.vector(x %*% b0)
        r_eta <- a * sqrt(2 * ub_obj)
        lb_eta <- y - r_eta
        ub_eta <- y + r_eta
        m_l_eta <- sum(abs(y)) + a * sqrt(2 * n * ub_obj)

        model$obj <- c(-xty, rep(0, p + 3 * n))
        model$Q <- Matrix::sparseMatrix(
            i = (2 * p + 1):(2 * p + n), j = (2 * p + 1):(2 * p + n),
            x = 0.5, dims = c(nv, nv)
        )

        # Rows: [1] sum(z) >= p - k;
        # [2..n+1] X beta - eta = 0;
        # [n+2..2n+1] eta - eta+ + eta- = 0;
        # [2n+2] sum(eta+ + eta-) <= M_l_eta.
        model$A <- Matrix::sparseMatrix(
            i = c(rep(1, p),
                  rep(2:(n + 1), p), 2:(n + 1),
                  rep((n + 2):(2 * n + 1), 3),
                  rep(2 * n + 2, 2 * n)),
            j = c((p + 1):(2 * p),
                  rep(1:p, each = n), (2 * p + 1):(2 * p + n),
                  (2 * p + 1):(2 * p + n),
                  (2 * p + n + 1):(2 * p + 2 * n),
                  (2 * p + 2 * n + 1):(2 * p + 3 * n),
                  (2 * p + n + 1):(2 * p + 3 * n)),
            x = c(rep(1, p),
                  as.vector(x), rep(-1, n),
                  rep(1, n), rep(-1, n), rep(1, n),
                  rep(1, 2 * n)),
            dims = c(2 * n + 2, nv)
        )
        model$rhs <- c(p - k, rep(0, 2 * n), m_l_eta)
        model$sense <- c(">", rep("=", 2 * n), "<")

        model$lb <- c(rep(-Inf, p), rep(0, p), lb_eta, rep(0, 2 * n))
        model$ub <- c(rep(Inf, p), rep(1, p), ub_eta,
                      pmax(ub_eta, 0), pmax(-lb_eta, 0))
        model$vtype <- c(rep("C", p), rep("B", p), rep("C", 3 * n))
        model$start <- c(b0, z0, eta0, pmax(eta0, 0), pmax(-eta0, 0))
    }

    # SOS-1 on (beta_i, z_i): at most one nonzero, so z_i = 1 =>
    # beta_i = 0.
    model$sos <- lapply(seq_len(p), function(i) {
        list(type = 1, index = c(i, p + i), weight = c(1, 2))
    })

    # MIPGap tightened from Gurobi's default 1e-4: nearly-tied subsets
    # differ by less than the default relative gap, and these MIOs are
    # small enough that full certification is cheap.
    params <- list(OutputFlag = as.numeric(verbose),
                   TimeLimit = time_limit, MIPGap = 1e-8)
    res <- gurobi::gurobi(model, params = params)

    beta <- res$x[seq_len(p)]
    # SOS-deselected coefficients are exact zeros up to solver tolerance;
    # impose them so that support counts are clean.
    beta[res$x[(p + 1):(2 * p)] > 0.5] <- 0

    list(beta = beta, status = res$status,
         mipgap = if (is.null(res$mipgap)) NA_real_ else res$mipgap,
         objval = res$objval)
}


# Information criteria on the residual mean squared error, matching the
# conventions used in train_lasso().
ic_eval <- function(mse, df, n, ic_type = "bic") {
    switch(ic_type,
        "aic" = n * log(mse) + 2 * df,
        "bic" = n * log(mse) + log(n) * df,
        "aicc" = n * log(mse) + 2 * df + 2 * df * (df + 1) / (n - df - 1),
        "hqc" = n * log(mse) + 2 * df * log(log(n)),
        stop("Unknown ic_type: ", ic_type)
    )
}
