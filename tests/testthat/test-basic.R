test_that("roll_predict works with basic methods", {
    n <- 50
    p <- 5
    x <- matrix(rnorm(n * p), n, p)
    colnames(x) <- paste0("x", 1:p)
    y <- x %*% rep(1, p) + rnorm(n)

    # Basic methods
    methods <- c("RW", "OLS", "PLasso")

    res <- roll_predict(x, y, roll_window = 30, h = 1, methods_use = methods, verbose = FALSE)

    expect_type(res, "list")
    expect_equal(names(res$mse), methods)
    expect_true(all(res$mse >= 0))

    # Check OLS dimensions
    expect_equal(length(res$OLS$y_hat), n - 30)
})

test_that("lasso returns correct structure", {
    x <- matrix(rnorm(50 * 5), 50, 5)
    y <- rnorm(50)

    mod <- lasso(x, y, method = "lasso")
    expect_s3_class(mod, "lasforecast_model")
    expect_true("lambda" %in% names(mod))
    expect_true("coefficients" %in% names(mod))
})
