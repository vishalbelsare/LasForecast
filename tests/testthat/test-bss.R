# Best subset selection via MIO (Bertsimas, King & Mazumder, 2016).
# Tests that touch the Gurobi solver are skipped when the (non-CRAN)
# 'gurobi' package is unavailable. Problem sizes are kept small so each
# MIO solves in well under a second.

best_subset_rss <- function(y, x, k, intercept = TRUE) {
    # Gold standard: exhaustive enumeration of all size-k subsets.
    combos <- utils::combn(ncol(x), k, simplify = FALSE)
    rss <- vapply(combos, function(S) {
        sum(stats::lsfit(x[, S, drop = FALSE], y, intercept = intercept)$residuals^2)
    }, numeric(1))
    list(rss = min(rss), support = combos[[which.min(rss)]])
}

test_that("first-order algorithm returns a k-sparse coefficient vector", {
    set.seed(1)
    n <- 40
    p <- 8
    x <- matrix(rnorm(n * p), n, p)
    y <- 2 * x[, 2] - x[, 5] + rnorm(n)
    xc <- sweep(x, 2, colMeans(x))
    yc <- y - mean(y)

    for (k in c(1, 3)) {
        b <- LasForecast:::bss_first_order(yc, xc, k)
        expect_length(b, p)
        expect_lte(sum(b != 0), k)
        # No worse than the all-zero solution
        expect_lte(sum((yc - xc %*% b)^2), sum(yc^2))
    }

    # High-dimensional initialization branch (p >= n)
    p2 <- 25
    n2 <- 20
    x2 <- matrix(rnorm(n2 * p2), n2, p2)
    y2 <- x2[, 1] + rnorm(n2)
    b2 <- LasForecast:::bss_first_order(y2, x2, 2)
    expect_lte(sum(b2 != 0), 2)
})

test_that("bss MIO matches exhaustive search (p < n, formulation 2.5)", {
    skip_if_not_installed("gurobi")
    set.seed(42)
    n <- 50
    p <- 6
    x <- matrix(rnorm(n * p), n, p)
    y <- 1 + 2 * x[, 1] - 1.5 * x[, 4] + rnorm(n)

    for (k in 1:2) {
        mod <- bss(x, y, k = k, time_limit = 30)
        cf <- coef(mod)
        expect_length(cf, p + 1)
        expect_lte(sum(cf[-1] != 0), k)

        rss_mio <- sum((y - cbind(1, x) %*% cf)^2)
        ex <- best_subset_rss(y, x, k)
        expect_lte(rss_mio, ex$rss + 1e-6)
        expect_identical(which(cf[-1] != 0), as.integer(ex$support))
    }
})

test_that("bss MIO matches exhaustive search (p >= n, formulation 2.6)", {
    skip_if_not_installed("gurobi")
    set.seed(7)
    n <- 25
    p <- 30
    x <- matrix(rnorm(n * p), n, p)
    y <- 0.5 + 2 * x[, 3] - 2 * x[, 17] + 0.5 * rnorm(n)

    mod <- bss(x, y, k = 2, time_limit = 30)
    cf <- coef(mod)
    rss_mio <- sum((y - cbind(1, x) %*% cf)^2)
    ex <- best_subset_rss(y, x, 2)
    expect_lte(rss_mio, ex$rss + 1e-6)
    expect_identical(which(cf[-1] != 0), as.integer(ex$support))
})

test_that("bss selects k by IC and returns a standard lasforecast_model", {
    skip_if_not_installed("gurobi")
    set.seed(3)
    n <- 50
    p <- 5
    x <- matrix(rnorm(n * p), n, p)
    y <- 1 + 2 * x[, 1] - 2 * x[, 3] + 0.5 * rnorm(n)

    mod <- bss(x, y, time_limit = 30)
    expect_s3_class(mod, "lasforecast_model")
    expect_identical(mod$method, "bss")
    expect_length(mod$coefficients, p + 1)

    # IC search path covers k = 0 (include_k0 default) through k_max = p
    expect_identical(mod$fit$ic_path$k, 0:p)
    expect_identical(mod$fit$k, mod$tuning_param)
    expect_identical(mod$fit$ic_path$ic[mod$fit$ic_path$k == mod$fit$k],
                     min(mod$fit$ic_path$ic))

    # Signal is strong: BIC should find the true support
    expect_identical(which(mod$coefficients[-1] != 0), c(1L, 3L))

    # S3 methods shared with lasso() outputs
    expect_identical(coef(mod), mod$coefficients)
    pred <- predict(mod, newx = x[1:3, ])
    expect_length(pred, 3)
    expect_output(print(mod), "Method: bss")

    # k = 0 excluded on request
    mod0 <- bss(x, y, include_k0 = FALSE, time_limit = 30)
    expect_identical(mod0$fit$ic_path$k, 1:p)
})

test_that("bss is robust to heterogeneous predictor scales", {
    skip_if_not_installed("gurobi")
    # Audit regression case: with columns on scales 100 vs 0.01, a
    # heuristic coefficient box computed from the warm start can cut off
    # the global optimum. The MIO must stay exact on the raw scale (no
    # rescaling of the data): only provable level-set bounds are imposed.
    set.seed(3)
    n <- 45
    p <- 5
    x <- matrix(rnorm(n * p), n, p)
    x <- sweep(x, 2, c(100, rep(1, p - 2), 0.01), "*")
    y <- as.numeric(0.3 + 0.02 * x[, 1] + 50 * x[, p] + rnorm(n))

    for (k in 2:3) {
        mod <- bss(x, y, k = k, time_limit = 30)
        cf <- coef(mod)
        rss <- sum((y - cbind(1, x) %*% cf)^2)
        expect_lte(rss, best_subset_rss(y, x, k)$rss + 1e-6)
    }
})

test_that("bss accepts k = 0 and a vector of candidate k values", {
    skip_if_not_installed("gurobi")
    set.seed(40)
    n <- 50
    p <- 5
    x <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(1 + 2 * x[, 1] - 2 * x[, 3] + 0.5 * rnorm(n))

    # Scalar k = 0: the null (prevailing-mean) model, no solver call
    m0 <- bss(x, y, k = 0)
    expect_identical(unname(coef(m0)), c(mean(y), rep(0, p)))
    expect_identical(m0$tuning_param, 0L)
    expect_null(m0$fit$ic_path)

    # Vector candidate set: IC selection restricted to the set
    mv <- bss(x, y, k = c(0, 2, 4), time_limit = 30)
    expect_identical(mv$fit$ic_path$k, c(0L, 2L, 4L))
    expect_true(mv$tuning_param %in% c(0L, 2L, 4L))
    expect_identical(mv$fit$ic_path$ic[mv$fit$ic_path$k == mv$tuning_param],
                     min(mv$fit$ic_path$ic))
    expect_lte(sum(coef(mv)[-1] != 0), mv$tuning_param)

    # Out-of-range candidates rejected
    expect_error(bss(x, y, k = c(1, p + 1)), "must be NULL")
    expect_error(bss(x, y, k = -1), "must be NULL")
    expect_error(bss(x, y, k = 1.5), "must be NULL")
})

test_that("bound multiplier a >= 1 leaves the exact solution unchanged", {
    skip_if_not_installed("gurobi")
    # a scales the provable level-set bound radii; any a >= 1 still contains
    # every candidate optimum, so the solution must be identical.
    set.seed(41)
    n <- 45
    p <- 6
    x <- sweep(matrix(rnorm(n * p), n, p), 2, c(100, rep(1, p - 2), 0.01), "*")
    y <- as.numeric(0.3 + 0.02 * x[, 1] + 50 * x[, p] + rnorm(n))

    m1 <- bss(x, y, k = 2, time_limit = 30)
    m2 <- bss(x, y, k = 2, a = 3, time_limit = 30)
    expect_equal(coef(m2), coef(m1), tolerance = 1e-6)

    # high-dimensional branch
    n2 <- 20
    p2 <- 25
    x2 <- matrix(rnorm(n2 * p2), n2, p2)
    y2 <- as.numeric(2 * x2[, 3] - 2 * x2[, 7] + 0.5 * rnorm(n2))
    h1 <- bss(x2, y2, k = 2, time_limit = 30)
    h2 <- bss(x2, y2, k = 2, a = 2, time_limit = 30)
    expect_equal(coef(h2), coef(h1), tolerance = 1e-6)

    expect_error(bss(x, y, k = 2, a = 0), "positive scalar")
})

test_that("bss is scale-equivariant: bss(xD, y) = D^-1 bss(x, y)", {
    skip_if_not_installed("gurobi")
    # The estimator is computed on the raw data with no internal
    # standardization, so rescaling columns must rescale coefficients
    # exactly inversely (the scale-equivariance of least squares on the
    # selected subset), leaving selection and fit unchanged.
    set.seed(8001)
    n <- 50
    p <- 6
    x <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(1 + 2 * x[, 1] - x[, 4] + rnorm(n))
    d <- 10^runif(p, -3, 3)

    for (k in c(1, 2)) {
        m1 <- bss(x, y, k = k, time_limit = 30)
        m2 <- bss(sweep(x, 2, d, "*"), y, k = k, time_limit = 30)
        cf1 <- coef(m1)
        cf2_back <- c(coef(m2)[1], coef(m2)[-1] * d)
        expect_equal(cf2_back, cf1, tolerance = 1e-6)
    }
})

test_that("constant predictor columns never receive a coefficient", {
    skip_if_not_installed("gurobi")
    # Rolling-window hazard: a predictor constant within the estimation
    # window is unidentified there but varies out of window, so a junk
    # coefficient would corrupt forecasts.
    set.seed(34)
    n <- 40
    p <- 4
    x <- matrix(rnorm(n * p), n, p)
    x[, 2] <- 1.7
    y <- as.numeric(0.5 + 2 * x[, 1] + rnorm(n))

    mod_ic <- bss(x, y, time_limit = 30)
    expect_identical(coef(mod_ic)[3], 0)
    mod_k <- bss(x, y, k = 2, time_limit = 30)
    expect_identical(coef(mod_k)[3], 0)
})

test_that("predictive = TRUE lag-aligns the data", {
    skip_if_not_installed("gurobi")
    set.seed(35)
    n <- 40
    x <- matrix(rnorm(n * 3), n, 3)
    y <- as.numeric(0.5 * x[, 1] + rnorm(n))

    mod_p <- bss(x, y, k = 1, predictive = TRUE, time_limit = 30)
    mod_m <- bss(x[-n, ], y[-1], k = 1, time_limit = 30)
    expect_identical(coef(mod_p), coef(mod_m))
})

test_that("bss handles intercept = FALSE and the k = p shortcut", {
    skip_if_not_installed("gurobi")
    set.seed(11)
    n <- 40
    p <- 4
    x <- matrix(rnorm(n * p), n, p)
    y <- x[, 2] + 0.5 * rnorm(n)

    mod <- bss(x, y, k = 1, intercept = FALSE, time_limit = 30)
    cf <- coef(mod)
    expect_length(cf, p)
    expect_lte(sum(cf != 0), 1)
    pred <- predict(mod, newx = x[1:2, ])
    expect_length(pred, 2)

    # k = p with p < n is unrestricted least squares
    mod_full <- bss(x, y, k = p, time_limit = 30)
    cf_ols <- unname(stats::lsfit(x, y)$coefficients)
    expect_equal(unname(coef(mod_full)), cf_ols, tolerance = 1e-6)
})

test_that("BSS integrates with roll_predict and backtest", {
    skip_if_not_installed("gurobi")
    set.seed(36)
    n <- 45
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- as.numeric(0.8 * x[, 1] + 0.3 * rnorm(n))

    # Fixed subset size per window
    res <- roll_predict(x, y, roll_window = 35, methods_use = c("RW", "BSS"),
                        bss_k = 1, verbose = FALSE, time_limit = 10)
    expect_length(res$BSS$y_hat, n - 35)
    expect_true(all(is.finite(res$BSS$y_hat)))
    expect_true(all(res$BSS$df <= 1))
    expect_true(all(res$BSS$tuning_param == 1))

    # IC-selected subset size: train_method = "bic" drives bss's ic_type
    bt <- backtest(x, y, methods = c("RWwD", "BSS"), roll_window = 35,
                   train_method = "bic", verbose = FALSE, time_limit = 10)
    expect_s3_class(bt, "lasforecast_backtest")
    rmse <- bt$summary_table$RMSE
    names(rmse) <- bt$summary_table$Method
    expect_true(is.finite(rmse[["BSS"]]))
    # The signal is strong, so BSS must beat the prevailing-mean benchmark
    expect_lt(rmse[["BSS"]], rmse[["RWwD"]])
})
