test_that("backtest works with simulated data", {
    set.seed(42)
    n <- 50
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- x %*% c(1, -0.5, 0) + rnorm(n)

    race <- backtest(x, y,
        methods = c("RW", "OLS", "PLasso"),
        roll_window = 30, h = 1, verbose = FALSE
    )

    expect_s3_class(race, "lasforecast_backtest")
    expect_true("summary_table" %in% names(race))
    expect_true("results" %in% names(race))
    expect_true("RMSE" %in% names(race$summary_table))
    expect_equal(nrow(race$summary_table), 4) # 3 methods + RWwD benchmark
})

test_that("backtest works with shift_y for multi-step forecasting", {
    set.seed(42)
    n <- 50
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- x %*% c(1, 0, 0) + rnorm(n)

    race <- backtest(x, y,
        methods = c("RW", "OLS"),
        roll_window = 25, h = 3, shift_y = TRUE, verbose = FALSE
    )

    expect_s3_class(race, "lasforecast_backtest")
    expect_equal(nrow(race$summary_table), 3) # 2 methods + RWwD benchmark
    expect_s3_class(race$results, "lasforecast_roll")
})

test_that("backtest works with data.frame input", {
    set.seed(42)
    df <- data.frame(
        x1 = rnorm(40),
        x2 = rnorm(40),
        target = rnorm(40)
    )

    race <- backtest(df, y = "target",
        methods = c("RW", "OLS"),
        roll_window = 20, verbose = FALSE
    )

    expect_s3_class(race, "lasforecast_backtest")
})

test_that("calc_loss works for all loss types", {
    actual <- c(1, 2, 3, 4, 5)
    predicted <- c(1.1, 2.2, 2.8, 4.1, 5.5)

    expect_true(is.numeric(calc_loss(actual, predicted, "mse")))
    expect_true(is.numeric(calc_loss(actual, predicted, "rmse")))
    expect_true(is.numeric(calc_loss(actual, predicted, "mae")))
    expect_true(is.numeric(calc_loss(actual, predicted, "mape")))
    expect_true(is.numeric(calc_loss(actual, predicted, "r2oos",
                                     benchmark_pred = rep(mean(actual), length(actual)))))
    expect_true(calc_loss(actual, predicted, "rmse") ==
                sqrt(calc_loss(actual, predicted, "mse")))
})

test_that("r2oos is the Campbell-Thompson out-of-sample R^2", {
    actual    <- c(1, 2, 3, 4, 5)
    predicted <- c(1.1, 2.2, 2.8, 4.1, 5.5)
    benchmark <- c(2, 2, 3, 3, 4)

    # No benchmark forecast supplied: the OOS R^2 is undefined -> NA
    expect_true(is.na(calc_loss(actual, predicted, "r2oos")))

    # Equals 1 - SSE_method / SSE_benchmark
    expect_equal(
        calc_loss(actual, predicted, "r2oos", benchmark_pred = benchmark),
        1 - sum((actual - predicted)^2) / sum((actual - benchmark)^2)
    )

    # The benchmark scored against itself has R2_OOS = 0
    expect_equal(
        calc_loss(actual, benchmark, "r2oos", benchmark_pred = benchmark),
        0
    )
})

test_that("calc_loss drops NA forecasts (failed windows)", {
    actual    <- c(1, 2, 3, 4, 5)
    predicted <- c(1.1, NA, 2.8, 4.1, 5.5)   # window 2 failed -> NA
    benchmark <- c(2, 2, 3, 3, 4)

    # Absolute metrics simply exclude the NA window
    expect_equal(calc_loss(actual, predicted, "mse"),
                 mean((actual - predicted)^2, na.rm = TRUE))
    expect_false(is.na(calc_loss(actual, predicted, "mae")))

    # r2oos sums numerator and denominator over the SAME (non-NA) windows
    ok <- !is.na(predicted)
    expect_equal(
        calc_loss(actual, predicted, "r2oos", benchmark_pred = benchmark),
        1 - sum((actual - predicted)[ok]^2) / sum((actual - benchmark)[ok]^2)
    )
})

test_that("alasso_model choice affects ALasso and TALasso results", {
    set.seed(123)
    n <- 60
    p <- 5
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    # Genuinely predictive target: y[t] is driven by x[t-1], matching the rolling
    # loop's y[t+1] ~ x[t] pairing, so the methods have real signal to select.
    # (With a contemporaneous y = x'beta + e there is no predictive signal and
    # BIC correctly returns the empty model in every window.)
    y <- c(0, as.numeric(x[-n, ] %*% c(2, -1, 0.5, 0, 0))) + rnorm(n)

    race_lasso <- backtest(x, y,
        methods = c("RW", "ALasso", "TALasso"),
        roll_window = 35, h = 1, verbose = FALSE,
        train_method = "bic", alasso_model = "lasso"
    )

    race_ridge <- backtest(x, y,
        methods = c("RW", "ALasso", "TALasso"),
        roll_window = 35, h = 1, verbose = FALSE,
        train_method = "bic", alasso_model = "ridge"
    )

    race_ols <- backtest(x, y,
        methods = c("RW", "ALasso", "TALasso"),
        roll_window = 35, h = 1, verbose = FALSE,
        train_method = "bic", alasso_model = "ols"
    )

    res_lasso <- race_lasso$results
    res_ridge <- race_ridge$results
    res_ols   <- race_ols$results

    # Different initial models should produce different beta_hat
    expect_false(isTRUE(all.equal(res_lasso$ALasso$beta_hat, res_ridge$ALasso$beta_hat)))
    expect_false(isTRUE(all.equal(res_lasso$ALasso$beta_hat, res_ols$ALasso$beta_hat)))

    # TALasso refines ALasso in its second adaptive step. When both recover the
    # same support their nonzero counts coincide, so the meaningful check is that
    # the coefficients differ, not the df.
    expect_false(isTRUE(all.equal(res_lasso$ALasso$beta_hat, res_lasso$TALasso$beta_hat)))
})

test_that("alasso_lambda_init flows through backtest to every window", {
    set.seed(5)
    n <- 100
    p <- 6
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- as.numeric(x %*% c(1.5, -1, rep(0, p - 2)) + rnorm(n))

    bt_fix <- backtest(x, y, methods = "ALasso", roll_window = 60, h = 1,
                       verbose = FALSE, train_method = "bic",
                       alasso_model = "lasso", alasso_lambda_init = 1e6)
    bt_def <- backtest(x, y, methods = "ALasso", roll_window = 60, h = 1,
                       verbose = FALSE, train_method = "bic",
                       alasso_model = "lasso")

    # The extreme initial penalty empties the model in every window.
    expect_true(all(bt_fix$results$ALasso$beta_hat == 0))
    expect_false(identical(bt_fix$results$ALasso$beta_hat,
                           bt_def$results$ALasso$beta_hat))
})

test_that("post_PLasso uses the raw scale and post variants honor alasso_model", {
    set.seed(7)
    n <- 70
    p <- 15
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    # Heterogeneous predictor scales + a genuinely predictive target (y[t] driven
    # by x[t-1]). Two effects need to show through: (i) raw-scale PLasso and
    # standardized SLasso select different supports, and (ii) with p moderately
    # large the OLS init overfits, so ols/ridge initial estimators yield
    # different adaptive supports. Both make the post-OLS refits differ. On
    # unit-scale, low-dim, or signal-free data the methods correctly coincide
    # (same support, or the empty model): that is not a scale bug.
    x[, 1] <- x[, 1] * 20
    x[, 2] <- x[, 2] * 0.1
    y <- c(0, as.numeric(x[-n, ] %*% c(0.1, -15, 1, rep(0, p - 3)))) + rnorm(n)

    # post_PLasso used to silently inherit scale_x = TRUE and coincide with
    # post_SLasso; with the raw-scale pin the two must differ.
    bt <- backtest(x, y, methods = c("post_PLasso", "post_SLasso"),
                   roll_window = 40, h = 1, verbose = FALSE, train_method = "bic")
    expect_false(identical(bt$results$post_PLasso$beta_hat,
                           bt$results$post_SLasso$beta_hat))

    # alasso_model now reaches the post variants (it changes the selected
    # support, which the post-OLS then refits).
    bt_ols <- backtest(x, y, methods = "post_ALasso", roll_window = 40, h = 1,
                       verbose = FALSE, train_method = "bic", alasso_model = "ols")
    bt_rdg <- backtest(x, y, methods = "post_ALasso", roll_window = 40, h = 1,
                       verbose = FALSE, train_method = "bic", alasso_model = "ridge")
    expect_false(identical(bt_ols$results$post_ALasso$beta_hat,
                           bt_rdg$results$post_ALasso$beta_hat))
})

test_that("ar_order = 1 forecasts are aligned and complete", {
    # The slot index used to be off by one for ar_order = 1: the first
    # feasible forecast was written to index 0 (a silent no-op), every stored
    # value shifted by one origin, and the last slot kept its initialization.
    set.seed(99)
    n <- 90
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- as.numeric(x %*% c(1.5, -1, 0) + rnorm(n))
    rw <- 60

    r <- suppressWarnings(roll_predict(x, y, roll_window = rw, h = 1,
        methods_use = "OLS", ar_order = 1, verbose = FALSE))
    nf <- length(r$OLS$y_hat)
    expect_equal(nf, n - rw - 1)

    # First feasible origin (i = rw + 2) must occupy slot 1...
    i <- rw + 2
    cf <- lsfit(cbind(x[(i - rw - 1):(i - 2), ], y[(i - rw - 1):(i - 2)]),
                y[(i - rw):(i - 1)])$coef
    expect_equal(r$OLS$y_hat[1], as.numeric(sum(c(1, c(x[i - 1, ], y[i - 1])) * cf)))

    # ...and the last origin (i = n) must fill the final slot.
    i <- n
    cf <- lsfit(cbind(x[(i - rw - 1):(i - 2), ], y[(i - rw - 1):(i - 2)]),
                y[(i - rw):(i - 1)])$coef
    expect_equal(r$OLS$y_hat[nf], as.numeric(sum(c(1, c(x[i - 1, ], y[i - 1])) * cf)))
})

test_that("expanding window with AR term supports h > 1", {
    # The expanding AR(1) block paired training rows starting at t = 2 with an
    # h-lag of y, which only exists from t = h + 1: cbind() crashed for h > 1.
    set.seed(99)
    n <- 80
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- as.numeric(x %*% c(1.5, -1, 0) + rnorm(n))

    r <- suppressWarnings(roll_predict(x, y, roll_window = 50, h = 2,
        methods_use = "OLS", ar_order = 1, window_type = "expanding",
        verbose = FALSE))
    expect_true(all(is.finite(r$OLS$y_hat)))
    expect_true(all(r$OLS$y_hat != 0))
})

test_that("wrappers reject reserved arguments instead of failing silently", {
    set.seed(1)
    n <- 60
    p <- 4
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- as.numeric(x %*% c(2, -1, 0, 0) + rnorm(n))

    # These used to error inside the rolling loop, get demoted to a warning,
    # and silently return all-NA forecasts for the method.
    expect_error(
        backtest(x, y, methods = "ALasso", roll_window = 40, verbose = FALSE,
                 train_method = "bic", model = "ols"),
        "alasso_model"
    )
    expect_error(
        backtest(x, y, methods = "ALasso", roll_window = 40, verbose = FALSE,
                 train_method = "bic", lambda_init = 0.1),
        "alasso_lambda_init"
    )
    expect_error(
        backtest(x, y, methods = "ALasso", roll_window = 40, verbose = FALSE,
                 train_method = "bic", predictive = TRUE),
        "predictive"
    )
    expect_error(
        backtest(x, y, methods = c("SLasso", "ALasso"), roll_window = 40,
                 verbose = FALSE, train_method = "bic", scale_x = TRUE),
        "scale_x"
    )

    # scale_x remains legitimate when no scale-pinned method is requested.
    bt <- backtest(x, y, methods = "ALasso", roll_window = 40, verbose = FALSE,
                   train_method = "bic", scale_x = FALSE)
    expect_s3_class(bt, "lasforecast_backtest")
})

test_that("print and summary methods work for lasforecast_backtest", {
    set.seed(42)
    n <- 40
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- rnorm(n)

    race <- backtest(x, y,
        methods = c("RW", "OLS"),
        roll_window = 20, verbose = FALSE
    )

    expect_output(print(race), "Backtest")

    summ <- summary(race)
    expect_true("RMSE_Ratio" %in% names(summ))
})
