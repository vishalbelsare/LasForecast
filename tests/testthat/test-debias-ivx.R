test_that("xdlasso returns correct class", {
    set.seed(123)
    n <- 100
    p <- 10
    w <- matrix(rnorm(n * p), n, p)
    colnames(w) <- paste0("x", 1:p)
    y <- w %*% c(0.5, -0.3, rep(0, p - 2)) + rnorm(n)

    fit <- xdlasso(w, y, d_ind = c(1, 2))

    expect_s3_class(fit, "lasforecast_debias_ivx")
    expect_equal(fit$d_ind, c(1, 2))
    expect_equal(fit$term_names, colnames(w))
    expect_length(fit$theta_hat_ivx, 2)
    expect_length(fit$sigma_hat_ivx, 2)
})

test_that("coef returns named focal estimates", {
    set.seed(123)
    n <- 100
    p <- 5
    w <- matrix(rnorm(n * p), n, p)
    colnames(w) <- paste0("x", 1:p)
    y <- w %*% c(1, 0, 0, 0, 0) + rnorm(n)

    fit <- xdlasso(w, y, d_ind = 1)
    cc <- coef(fit)

    expect_length(cc, 1)
    expect_equal(names(cc), "x1")
})

test_that("summary produces coefficient table with p-values", {
    set.seed(42)
    n <- 100
    p <- 8
    w <- matrix(rnorm(n * p), n, p)
    colnames(w) <- paste0("x", 1:p)
    y <- w %*% c(0.8, 0, -0.5, rep(0, p - 3)) + rnorm(n)

    fit <- xdlasso(w, y, d_ind = c(1, 3))
    s <- summary(fit)

    expect_s3_class(s, "summary.lasforecast_debias_ivx")
    expect_equal(nrow(s$focal), 2)
    expect_equal(colnames(s$focal), c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    expect_true(all(s$focal[, "Std. Error"] > 0))
    expect_true(all(s$focal[, "Pr(>|t|)"] >= 0 & s$focal[, "Pr(>|t|)"] <= 1))
})

test_that("summary includes lasso controls", {
    set.seed(42)
    n <- 100
    p <- 8
    w <- matrix(rnorm(n * p), n, p)
    colnames(w) <- paste0("x", 1:p)
    y <- w %*% c(0.8, 0.6, -0.5, rep(0, p - 3)) + rnorm(n)

    fit <- xdlasso(w, y, d_ind = 1)
    s <- summary(fit)

    # x2 and x3 should likely appear as lasso-selected controls
    expect_true(is.numeric(s$lasso_controls))
})

test_that("print methods run without error", {
    set.seed(99)
    n <- 80
    p <- 5
    w <- matrix(rnorm(n * p), n, p)
    colnames(w) <- paste0("x", 1:p)
    y <- w %*% c(1, 0, 0, 0, 0) + rnorm(n)

    fit <- xdlasso(w, y, d_ind = 1)

    expect_output(print(fit), "IVX-Desparsified")
    expect_output(print(summary(fit)), "Focal coefficients")
})

test_that("single-restriction Wald equals the squared t-statistic", {
    set.seed(42)
    n <- 150
    p <- 6
    w <- matrix(rnorm(n * p), n, p)
    colnames(w) <- paste0("x", 1:p)
    y <- as.numeric(w %*% c(0.5, -0.3, rep(0, p - 2)) + rnorm(n))

    # The joint covariance must follow the marginal-SE convention, so for a
    # single restriction the Wald stat reduces to (theta_hat / se)^2 under both
    # the iid (default) and robust standard errors.
    for (st in c("iid", "robust")) {
        fit <- debias_ivx(w, y, d_ind = 1, joint_test = TRUE, se_type = st)
        t2 <- (fit$theta_hat_ivx[1] / fit$sigma_hat_ivx[1])^2
        expect_equal(fit$wald_stat, t2)
    }
})

test_that("joint Wald test shown in summary", {
    set.seed(42)
    n <- 100
    p <- 5
    w <- matrix(rnorm(n * p), n, p)
    colnames(w) <- paste0("x", 1:p)
    y <- w %*% c(0.5, -0.3, rep(0, p - 2)) + rnorm(n)

    fit <- xdlasso(w, y, d_ind = c(1, 2), joint_test = TRUE)
    s <- summary(fit)

    expect_false(is.na(s$wald_stat))
    expect_false(is.na(s$p_value_wald))
    expect_output(print(s), "Wald test")
})

test_that("works with unnamed matrix", {
    set.seed(123)
    n <- 80
    p <- 5
    w <- matrix(rnorm(n * p), n, p)
    y <- w %*% c(1, rep(0, p - 1)) + rnorm(n)

    fit <- xdlasso(w, y, d_ind = 1)

    expect_equal(fit$term_names, paste0("X", 1:p))
    expect_equal(names(coef(fit)), "X1")
})

test_that("data.frame input with d (character names)", {
    set.seed(10)
    n <- 80
    p <- 10
    df <- data.frame(
        date = seq.Date(as.Date("2000-01-01"), by = "month", length.out = n),
        y = rnorm(n),
        matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
    )
    fit <- xdlasso(df, y = "y", d = "x1", train_method = "bic")
    expect_s3_class(fit, "lasforecast_debias_ivx")
    expect_equal(fit$term_names[fit$d_ind], "x1")
})

test_that("predictive auto-detects from date presence", {
    set.seed(20)
    n <- 80
    p <- 10
    pred <- matrix(rnorm(n * p), n, p, dimnames = list(NULL, paste0("x", 1:p)))
    y_vec <- rnorm(n)
    df <- data.frame(
        date = seq.Date(as.Date("2000-01-01"), by = "month", length.out = n),
        y = y_vec, pred
    )
    # data.frame with date: predictive=TRUE by default
    fit_df <- xdlasso(df, y = "y", d = "x1", train_method = "bic")

    # Manual pre-aligned call should match auto-aligned
    w2 <- pred[-n, , drop = FALSE]
    y2 <- y_vec[-1]
    fit_manual <- xdlasso(w2, y2, d_ind = 1, predictive = FALSE, train_method = "bic")
    expect_equal(fit_df$theta_hat_ivx, fit_manual$theta_hat_ivx)
})

test_that("d vs d_ind validation", {
    w <- matrix(rnorm(200), 40, 5)
    colnames(w) <- paste0("x", 1:5)
    y <- rnorm(40)

    expect_error(xdlasso(w, y, d = "x1", d_ind = 1), "not both")
    expect_error(xdlasso(w, y), "must be specified")
    expect_error(xdlasso(w, y, d = "MISSING"), "not found")
})
