test_that("arma() with a fixed order returns a valid lasforecast_model", {
    set.seed(1)
    y <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.3), n = 150))

    mod <- arma(y = y, order = c(1, 0, 1), h = 1)

    expect_s3_class(mod, "lasforecast_model")
    expect_equal(mod$method, "arma")
    expect_equal(mod$arma_order, c(1, 0, 1))
    expect_true(is.finite(mod$forecast))
    expect_length(mod$forecast, 1)
    # ARMA has no predictor coefficients
    expect_true(all(is.na(mod$coefficients)))
})

test_that("arma() auto-selects an order when none is given", {
    set.seed(2)
    y <- as.numeric(arima.sim(list(ar = c(0.6, -0.2)), n = 200))

    mod <- arma(y = y, max_p = 3, max_q = 3)

    expect_s3_class(mod, "lasforecast_model")
    expect_true(is.finite(mod$forecast))
    expect_length(mod$arma_order, 3)
    expect_true(mod$df >= 1)
})

test_that("predict() returns the precomputed forecast and ignores newx", {
    set.seed(3)
    y <- as.numeric(arima.sim(list(ar = 0.4), n = 120))
    mod <- arma(y = y, order = c(1, 0, 0), h = 1)

    # newx is irrelevant for a univariate model
    expect_equal(predict(mod, newx = c(1, 2, 3)), mod$forecast)
    expect_equal(predict(mod, newx = matrix(rnorm(10), 1)), mod$forecast)
})

test_that("h-step forecast aligns with stats::arima predict", {
    set.seed(4)
    y <- as.numeric(arima.sim(list(ar = 0.5, ma = 0.2), n = 150))

    mod <- arma(y = y, order = c(1, 0, 1), h = 3)
    ref <- stats::arima(y, order = c(1, 0, 1), include.mean = TRUE)
    expect_equal(mod$forecast, as.numeric(predict(ref, n.ahead = 3)$pred)[3])
})

test_that("backtest runs with ARMA as a benchmark", {
    set.seed(42)
    n <- 80
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- as.numeric(arima.sim(list(ar = 0.5), n = n))

    bt <- backtest(x, y,
        methods = c("RW", "OLS", "ARMA"),
        roll_window = 50, h = 1, verbose = FALSE
    )

    expect_s3_class(bt, "lasforecast_backtest")
    expect_true("ARMA" %in% bt$summary_table$Method)
    arma_rmse <- bt$summary_table$RMSE[bt$summary_table$Method == "ARMA"]
    expect_true(is.finite(arma_rmse) && arma_rmse > 0)
})

test_that("roll_predict ARMA produces real forecasts and NA betas", {
    set.seed(7)
    n <- 70
    x <- matrix(rnorm(n * 2), n, 2)
    colnames(x) <- c("x1", "x2")
    y <- as.numeric(arima.sim(list(ar = 0.4, ma = 0.3), n = n))

    res <- roll_predict(x, y,
        roll_window = 45, h = 1, methods_use = "ARMA",
        verbose = FALSE
    )

    # Forecasts are populated (not left at the initialized 0/NA)
    expect_true(all(is.finite(res$ARMA$y_hat)))
    # No predictor coefficients
    expect_true(all(is.na(res$ARMA$beta_hat)))
    expect_true(is.finite(res$mse["ARMA"]))
})

test_that("ARMA estimation specs flow through backtest's ...", {
    set.seed(11)
    n <- 80
    x <- matrix(rnorm(n * 2), n, 2)
    colnames(x) <- c("x1", "x2")
    y <- as.numeric(arima.sim(list(ar = 0.6), n = n))

    # Passing ic / max orders via ... must not error or collide with lasso args
    bt <- backtest(x, y,
        methods = c("RWwD", "ARMA"),
        roll_window = 50, verbose = FALSE,
        ic = "bic", max_p = 2, max_q = 2
    )
    expect_s3_class(bt, "lasforecast_backtest")
    expect_true(is.finite(bt$summary_table$RMSE[bt$summary_table$Method == "ARMA"]))
})

test_that("arma() falls back gracefully on a constant series", {
    y <- rep(3.5, 60)
    expect_no_error(mod <- arma(y = y, order = c(1, 0, 0)))
    expect_true(is.finite(mod$forecast))
})

test_that("print method works for an ARMA model", {
    set.seed(5)
    y <- as.numeric(arima.sim(list(ar = 0.5), n = 100))
    mod <- arma(y = y, order = c(1, 0, 0))
    expect_output(print(mod), "ARMA order")
})
