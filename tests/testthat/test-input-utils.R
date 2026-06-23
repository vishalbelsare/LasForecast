test_that("prepare_input works with matrix input", {
    x <- matrix(rnorm(50 * 5), 50, 5)
    colnames(x) <- paste0("x", 1:5)
    y <- rnorm(50)

    inp <- prepare_input(x, y)
    expect_true(is.matrix(inp$x_mat))
    expect_true(is.numeric(inp$y_vec))
    expect_equal(nrow(inp$x_mat), 50)
    expect_equal(ncol(inp$x_mat), 5)
    expect_equal(length(inp$y_vec), 50)
    expect_null(inp$dates)
    expect_equal(inp$var_names, paste0("x", 1:5))
})

test_that("prepare_input works with data.frame and y as column name", {
    df <- data.frame(
        date = seq(as.Date("2020-01-01"), by = "month", length.out = 50),
        x1 = rnorm(50),
        x2 = rnorm(50),
        target = rnorm(50)
    )

    inp <- prepare_input(df, y = "target")
    expect_equal(ncol(inp$x_mat), 2)
    expect_equal(length(inp$y_vec), 50)
    expect_false(is.null(inp$dates))
    expect_equal(length(inp$dates), 50)
    expect_true("x1" %in% inp$var_names)
    expect_false("target" %in% inp$var_names)
    expect_false("date" %in% inp$var_names)
})

test_that("prepare_input works with data.frame and y as vector", {
    df <- data.frame(x1 = rnorm(30), x2 = rnorm(30))
    y <- rnorm(30)

    inp <- prepare_input(df, y = y)
    expect_equal(ncol(inp$x_mat), 2)
    expect_equal(length(inp$y_vec), 30)
    expect_null(inp$dates)
})

test_that("prepare_input with explicit date_col", {
    df <- data.frame(
        my_date = seq(as.Date("2020-01-01"), by = "month", length.out = 20),
        x1 = rnorm(20),
        target = rnorm(20)
    )

    inp <- prepare_input(df, y = "target", date_col = "my_date")
    expect_false(is.null(inp$dates))
    expect_equal(ncol(inp$x_mat), 1)
})

test_that("validate_input catches dimension mismatch", {
    x <- matrix(rnorm(50 * 5), 50, 5)
    y <- rnorm(30)

    expect_error(prepare_input(x, y), "Dimension mismatch")
})

test_that("validate_input catches non-numeric", {
    x <- matrix(rnorm(20), 10, 2)
    y <- letters[1:10]

    expect_error(prepare_input(x, y), "numeric")
})

test_that("lasso works with data.frame input", {
    df <- data.frame(
        x1 = rnorm(50),
        x2 = rnorm(50),
        x3 = rnorm(50),
        target = rnorm(50)
    )

    mod <- lasso(df, y = "target", method = "lasso")
    expect_s3_class(mod, "lasforecast_model")
    expect_true("lambda" %in% names(mod))
})

test_that("roll_predict works with data.frame input", {
    df <- data.frame(
        x1 = rnorm(60),
        x2 = rnorm(60),
        target = rnorm(60)
    )

    res <- roll_predict(df, y = "target", roll_window = 30, h = 1,
                        methods_use = c("RW", "OLS"), verbose = FALSE)
    expect_s3_class(res, "lasforecast_roll")
    expect_equal(length(res$OLS$y_hat), 30)
})

test_that("lasso predictive auto-detects from date presence", {
    set.seed(30)
    n <- 80
    p <- 5
    pred <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
    y_vec <- rnorm(n)
    df <- data.frame(
        date = seq.Date(as.Date("2000-01-01"), by = "month", length.out = n),
        y = y_vec, pred
    )
    # data.frame with date: predictive=TRUE by default (n-1 effective obs)
    mod_auto <- lasso(df, y = "y", method = "lasso", train_method = "bic")

    # Manual lag alignment should match
    x_lag <- pred[-n, , drop = FALSE]
    y_lead <- y_vec[-1]
    mod_manual <- lasso(x_lag, y_lead, method = "lasso", train_method = "bic",
                        predictive = FALSE)
    expect_equal(mod_auto$coefficients, mod_manual$coefficients)
})

test_that("lasso predictive=FALSE with date data.frame skips alignment", {
    set.seed(31)
    n <- 80
    df <- data.frame(
        date = seq.Date(as.Date("2000-01-01"), by = "month", length.out = n),
        y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
    )
    # Explicit FALSE: no alignment, uses all n obs
    mod_no <- lasso(df, y = "y", method = "lasso", train_method = "bic",
                    predictive = FALSE)
    # Explicit TRUE: loses 1 obs
    mod_yes <- lasso(df, y = "y", method = "lasso", train_method = "bic",
                     predictive = TRUE)
    # Coefficients should differ (different training data)
    expect_false(isTRUE(all.equal(mod_no$coefficients, mod_yes$coefficients)))
})

test_that("roll_predict does not double-lag with predictive lasso", {
    set.seed(32)
    n <- 80
    df <- data.frame(
        date = seq.Date(as.Date("2000-01-01"), by = "month", length.out = n),
        y = rnorm(n), x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n)
    )
    # roll_predict passes raw matrices internally, so lasso predictive defaults FALSE
    res <- roll_predict(df, y = "y", roll_window = 40, h = 1,
                        methods_use = c("SLasso"), verbose = FALSE,
                        train_method = "bic")
    expect_s3_class(res, "lasforecast_roll")
    # Should produce n - roll_window = 40 forecasts (no extra obs lost)
    expect_equal(length(res$SLasso$y_hat), n - 40)
})
