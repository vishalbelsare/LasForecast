test_that("autoplot.lasforecast_roll returns ggplot", {
    set.seed(42)
    n <- 60
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- rnorm(n)

    res <- roll_predict(x, y,
        roll_window = 30, h = 1,
        methods_use = c("RW", "OLS", "PLasso"), verbose = FALSE
    )

    p1 <- autoplot(res)
    expect_s3_class(p1, "ggplot")
})

test_that("autoplot.lasforecast_roll uses dates when available", {
    set.seed(42)
    n <- 60
    p <- 3
    df <- data.frame(
        date = seq(as.Date("2020-01-01"), by = "month", length.out = n),
        x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
        target = rnorm(n)
    )

    res <- roll_predict(df, y = "target",
        roll_window = 30, h = 1,
        methods_use = c("RW", "OLS"), verbose = FALSE
    )

    expect_false(is.null(res$dates))
    p1 <- autoplot(res, type = "forecasts")
    expect_s3_class(p1, "ggplot")
})

test_that("autoplot.lasforecast_backtest returns ggplot for bar type", {
    set.seed(42)
    n <- 60
    p <- 3
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- rnorm(n)

    race <- backtest(x, y,
        methods = c("RW", "OLS", "PLasso"),
        roll_window = 30, verbose = FALSE
    )

    p1 <- autoplot(race, type = "bar")
    expect_s3_class(p1, "ggplot")
})

test_that("plot_variable_importance returns ggplot", {
    set.seed(42)
    n <- 60
    p <- 5
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- x %*% c(1, 0, 0, 0, 0) + rnorm(n)

    res <- roll_predict(x, y,
        roll_window = 30, h = 1,
        methods_use = c("RW", "PLasso"), verbose = FALSE
    )

    p1 <- plot_variable_importance(res, methods = "PLasso", type = "frequency")
    expect_s3_class(p1, "ggplot")

    p2 <- plot_variable_importance(res, methods = "PLasso", type = "heatmap")
    expect_s3_class(p2, "ggplot")

    p3 <- plot_variable_importance(res, methods = "PLasso", type = "evolution")
    expect_s3_class(p3, "ggplot")
})


test_that("compare_forecasts works with dates and configurable loss", {
    set.seed(42)
    n <- 60
    p <- 3
    df <- data.frame(
        date = seq(as.Date("2020-01-01"), by = "month", length.out = n),
        x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
        target = rnorm(n)
    )

    res <- roll_predict(df, y = "target",
        roll_window = 30, h = 1,
        methods_use = c("RW", "OLS"), verbose = FALSE
    )

    comp <- compare_forecasts(res, loss = c("rmse", "mae", "mape"))
    expect_true("MAPE" %in% names(comp))
    expect_true("RMSE" %in% names(comp))
})

test_that("theme_lasforecast returns a valid theme", {
    th <- theme_lasforecast()
    expect_s3_class(th, "theme")
})
