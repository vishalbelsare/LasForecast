#' Prepare input data for LasForecast functions
#'
#' Accepts matrix, data.frame, tibble, or tsibble inputs and converts them
#' to the standard format used internally (numeric matrix x, numeric vector y).
#'
#' @param x Predictor data. Can be a matrix, data.frame, tibble, or tsibble.
#'   When \code{x} is a data.frame/tibble/tsibble and \code{y} is a string,
#'   \code{y} is interpreted as a column name in \code{x}.
#' @param y Response vector, or a column name string when \code{x} is a data frame.
#'   Can be NULL if the response is included in \code{x} and specified by name.
#' @param date_col Character string specifying the date column name.
#'   For tsibble objects, the index column is used automatically.
#'   For data.frames/tibbles, the function auto-detects Date/POSIXct columns
#'   if \code{date_col} is NULL.
#'
#' @return A list with components:
#' \describe{
#'   \item{x_mat}{Numeric predictor matrix}
#'   \item{y_vec}{Numeric response vector}
#'   \item{dates}{Date vector (or NULL if no dates found)}
#'   \item{var_names}{Character vector of predictor variable names}
#' }
#'
#' @keywords internal
prepare_input <- function(x, y = NULL, date_col = NULL) {
    dates <- NULL
    var_names <- NULL

    # --- Handle tsibble ---
    if (inherits(x, "tbl_ts")) {
        # Extract index as dates
        idx_var <- tsibble::index_var(x)
        dates <- as.Date(x[[idx_var]])

        # Remove index and key columns
        key_vars <- tsibble::key_vars(x)
        drop_cols <- unique(c(idx_var, key_vars))

        # Convert to regular tibble first
        x_df <- as.data.frame(x)
        x_df <- x_df[, !names(x_df) %in% drop_cols, drop = FALSE]

        # Extract y if specified as string
        if (is.character(y) && length(y) == 1) {
            if (!y %in% names(x_df)) {
                stop("Column '", y, "' not found in the tsibble.")
            }
            y_vec <- as.numeric(x_df[[y]])
            x_df <- x_df[, !names(x_df) %in% y, drop = FALSE]
        } else if (is.null(y)) {
            stop("When x is a tsibble, y must be specified as a column name string.")
        } else {
            y_vec <- as.numeric(y)
        }

        var_names <- names(x_df)
        x_mat <- as.matrix(x_df)
        storage.mode(x_mat) <- "double"

    # --- Handle data.frame / tibble ---
    } else if (is.data.frame(x)) {
        x_df <- as.data.frame(x)

        # Auto-detect or use specified date column
        if (!is.null(date_col)) {
            if (!date_col %in% names(x_df)) {
                stop("date_col '", date_col, "' not found in x.")
            }
            dates <- as.Date(x_df[[date_col]])
            x_df <- x_df[, !names(x_df) %in% date_col, drop = FALSE]
        } else {
            # Auto-detect Date/POSIXct columns
            date_candidates <- vapply(x_df, function(col) {
                inherits(col, "Date") || inherits(col, "POSIXct") || inherits(col, "POSIXlt")
            }, logical(1))
            if (any(date_candidates)) {
                date_col_name <- names(x_df)[which(date_candidates)[1]]
                dates <- as.Date(x_df[[date_col_name]])
                x_df <- x_df[, !names(x_df) %in% date_col_name, drop = FALSE]
            }
        }

        # Extract y if specified as string
        if (is.character(y) && length(y) == 1) {
            if (!y %in% names(x_df)) {
                stop("Column '", y, "' not found in x.")
            }
            y_vec <- as.numeric(x_df[[y]])
            x_df <- x_df[, !names(x_df) %in% y, drop = FALSE]
        } else if (is.null(y)) {
            stop("y must be specified (as a vector or column name string).")
        } else {
            y_vec <- as.numeric(y)
        }

        var_names <- names(x_df)
        x_mat <- as.matrix(x_df)
        storage.mode(x_mat) <- "double"

    # --- Handle matrix ---
    } else {
        x_mat <- as.matrix(x)
        storage.mode(x_mat) <- "double"
        var_names <- colnames(x_mat)

        if (is.null(y)) {
            stop("y must be provided when x is a matrix.")
        }
        y_vec <- as.numeric(y)
    }

    # Assign default column names if missing
    if (is.null(var_names) || all(var_names == "")) {
        var_names <- paste0("X", seq_len(ncol(x_mat)))
        colnames(x_mat) <- var_names
    }

    # Validate
    validate_input(x_mat, y_vec)

    list(
        x_mat = x_mat,
        y_vec = y_vec,
        dates = dates,
        var_names = var_names
    )
}

#' Validate input dimensions and types
#'
#' @param x Numeric matrix
#' @param y Numeric vector
#'
#' @keywords internal
validate_input <- function(x, y) {
    if (!is.matrix(x) || !is.numeric(x)) {
        stop("x must be a numeric matrix after preparation.")
    }
    if (!is.numeric(y) || all(is.na(y))) {
        stop("y must be a numeric vector.")
    }
    if (nrow(x) != length(y)) {
        stop(
            "Dimension mismatch: x has ", nrow(x), " rows but y has ",
            length(y), " elements."
        )
    }
    if (anyNA(x)) {
        warning("x contains NA values. Consider removing or imputing them.")
    }
    if (anyNA(y)) {
        warning("y contains NA values. Consider removing or imputing them.")
    }
    invisible(TRUE)
}
