test_that("Advanced estimation methods work", {
    n <- 50
    p <- 5
    x <- matrix(rnorm(n * p), n, p)
    y <- x %*% c(1, rep(0, p - 1)) + rnorm(n)

    # ALasso
    alas <- lasso(x, y, method = "alasso")
    expect_s3_class(alas, "lasforecast_model")

    # OLS
    mod_ols <- ols(x, y)
    expect_s3_class(mod_ols, "lasforecast_model")

    # Post-Lasso
    plasso <- lasso(x, y, method = "post_lasso")
    expect_s3_class(plasso, "lasforecast_model")
})

test_that("ALasso with scale_x = TRUE is invariant to predictor rescaling", {
    # A standardized adaptive Lasso must not depend on the predictors' units:
    # its weights come from the initial estimate read on the standardized scale.
    set.seed(1)
    n <- 120
    p <- 8
    x <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(x %*% c(1.5, -1, rep(0, p - 2)) + rnorm(n))

    # Rescale columns to wildly different units (same information).
    s <- exp(seq(-3, 3, length.out = p))
    x2 <- sweep(x, 2, s, "*")

    f1 <- lasso(x, y, method = "alasso", predictive = FALSE,
                scale_x = TRUE, train_method = "cv")
    f2 <- lasso(x2, y, method = "alasso", predictive = FALSE,
                scale_x = TRUE, train_method = "cv")

    # Forecasts on the corresponding (rescaled) design must coincide.
    newx <- matrix(rnorm(5 * p), 5, p)
    p1 <- predict(f1, newx = newx)
    p2 <- predict(f2, newx = sweep(newx, 2, s, "*"))
    expect_equal(p1, p2, tolerance = 1e-6)
})

test_that("TALasso with scale_x = TRUE is invariant to predictor rescaling", {
    # Same requirement as for ALasso: the round-2 weights must be read on the
    # standardized scale when the round-2 penalty standardizes (talasso()).
    set.seed(1)
    n <- 120
    p <- 8
    x <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(x %*% c(1.5, -1, rep(0, p - 2)) + rnorm(n))

    s <- exp(seq(-3, 3, length.out = p))
    x2 <- sweep(x, 2, s, "*")

    f1 <- lasso(x, y, method = "talasso", predictive = FALSE,
                scale_x = TRUE, train_method = "cv")
    f2 <- lasso(x2, y, method = "talasso", predictive = FALSE,
                scale_x = TRUE, train_method = "cv")

    newx <- matrix(rnorm(5 * p), 5, p)
    p1 <- predict(f1, newx = newx)
    p2 <- predict(f2, newx = sweep(newx, 2, s, "*"))
    expect_equal(p1, p2, tolerance = 1e-6)
})

test_that("lasso_weight_single maps coefficients back to the original scale", {
    # With scale_x = TRUE the soft-threshold acts on the standardized slope;
    # the returned slope must be on the original scale (single-survivor path
    # of TALasso), so rescaling x by 1/10 must scale the slope by 10.
    set.seed(2)
    n <- 100
    x <- rnorm(n, sd = 5)
    y <- 0.8 * x + rnorm(n)

    r1 <- LasForecast:::lasso_weight_single(x, y, lambda = 0.05,
                                            w = 1, scale_x = TRUE)
    r2 <- LasForecast:::lasso_weight_single(x / 10, y, lambda = 0.05,
                                            w = 1, scale_x = TRUE)

    expect_true(r1$bhat != 0)
    expect_equal(r2$bhat, r1$bhat * 10)
    expect_equal(r2$ahat, r1$ahat)
})

test_that("TALasso warns when the first stage uses a Lasso initial estimator", {
    set.seed(3)
    n <- 80
    p <- 6
    x <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(x %*% c(1.5, -1, rep(0, p - 2)) + rnorm(n))

    expect_warning(
        lasso(x, y, method = "talasso", predictive = FALSE,
              train_method = "bic", model = "lasso"),
        "low-dimensional"
    )
})

test_that("Lasso-screened variables stay excluded under extreme units", {
    # The zero-floor in the adaptive weights must live on the penalty scale:
    # a variable screened out by the initial Lasso carries an effectively
    # infinite weight regardless of its units, so it cannot re-enter the
    # adaptive step just because its column is measured in huge units.
    set.seed(1)
    n <- 120
    p <- 25
    x <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(x %*% c(1.5, -1, rep(0, p - 2)) + rnorm(n))

    f1 <- lasso(x, y, method = "alasso", predictive = FALSE,
                train_method = "cv")
    zeroed <- which(f1$coefficients[-1] == 0)[1]
    s <- rep(1, p)
    s[zeroed] <- 1e12

    f2 <- lasso(sweep(x, 2, s, "*"), y, method = "alasso",
                predictive = FALSE, train_method = "cv")

    newx <- matrix(rnorm(5 * p), 5, p)
    p1 <- predict(f1, newx = newx)
    p2 <- predict(f2, newx = sweep(newx, 2, s, "*"))
    expect_equal(p1, p2, tolerance = 1e-6)
})

test_that("ALasso default scale_x follows the initial estimator", {
    set.seed(4)
    n <- 120
    xl <- matrix(rnorm(n * 4), n, 4)                   # p < 2*sqrt(n) -> OLS init
    yl <- as.numeric(xl %*% c(1.5, -1, 0, 0) + rnorm(n))
    xh <- matrix(rnorm(n * 25), n, 25)                 # p > 2*sqrt(n) -> Lasso init
    yh <- as.numeric(xh %*% c(1.5, -1, rep(0, 23)) + rnorm(n))

    # Low-dim default reproduces the explicit raw-scale fit (LSG 2022)...
    fl <- lasso(xl, yl, method = "alasso", predictive = FALSE,
                train_method = "bic")
    fl_raw <- lasso(xl, yl, method = "alasso", predictive = FALSE,
                    train_method = "bic", scale_x = FALSE)
    expect_equal(fl$coefficients, fl_raw$coefficients)

    # ...while the high-dim default keeps the standardized convention.
    fh <- lasso(xh, yh, method = "alasso", predictive = FALSE,
                train_method = "bic")
    fh_std <- lasso(xh, yh, method = "alasso", predictive = FALSE,
                    train_method = "bic", scale_x = TRUE)
    expect_equal(fh$coefficients, fh_std$coefficients)
})

test_that("lambda_init reaches both the tuning and the estimation steps", {
    set.seed(5)
    n <- 100
    p <- 6
    x <- matrix(rnorm(n * p), n, p)
    y <- as.numeric(x %*% c(1.5, -1, rep(0, p - 2)) + rnorm(n))

    # The public call must reproduce the manual composition with the same
    # initial-step penalty, so tuning and estimation share one set of weights.
    L <- 0.05
    f <- lasso(x, y, method = "alasso", predictive = FALSE, model = "lasso",
               scale_x = TRUE, train_method = "bic", lambda_init = L)
    lam <- train_lasso(x, y, ada = TRUE, model = "lasso", lambda_init = L,
                       intercept = TRUE, scale_x = TRUE, train_method = "bic")
    fit <- LasForecast:::adalasso(x, y, lambda = lam, lambda_init = L,
                                  model = "lasso", scale_x = TRUE)
    expect_identical(f$coefficients, c(fit$ahat, fit$bhat))

    # An extreme initial penalty screens out everything: the adaptive step
    # must return the empty model rather than error on a degenerate grid.
    f2 <- lasso(x, y, method = "alasso", predictive = FALSE, model = "lasso",
                scale_x = TRUE, train_method = "bic", lambda_init = 1e6)
    expect_true(all(f2$coefficients[-1] == 0))
})

test_that("Horse Racing works", {
    n <- 100
    p <- 5
    x <- matrix(rnorm(n * p), n, p)
    y <- x %*% c(1, rep(0, p - 1)) + rnorm(n)

    res <- roll_predict(x, y,
        roll_window = 50, h = 1,
        methods_use = c("OLS", "PLasso", "ALasso"), verbose = FALSE
    )

    # Compare
    comp <- compare_forecasts(res)
    expect_s3_class(comp, "data.frame")
    expect_true("RMSE" %in% names(comp))
})

test_that("Data is available", {
    expect_true(exists("fredmd"))
})
