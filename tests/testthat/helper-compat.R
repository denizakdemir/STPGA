if (!exists("expect_no_error", mode = "function")) {
  expect_no_error <- function(object, regexp = NA, ...) {
    testthat::expect_error(object, regexp = regexp, ...)
  }
}

expect_true <- function(object, info = NULL, label = NULL) {
  act <- substitute(object)
  value <- eval.parent(act)
  testthat::expect(
    isTRUE(value),
    sprintf("%s is not TRUE", if (is.null(label)) deparse(act) else label),
    info = info
  )
  invisible(value)
}

expect_false <- function(object, info = NULL, label = NULL) {
  act <- substitute(object)
  value <- eval.parent(act)
  testthat::expect(
    identical(value, FALSE),
    sprintf("%s is not FALSE", if (is.null(label)) deparse(act) else label),
    info = info
  )
  invisible(value)
}

expect_null <- function(object, info = NULL, label = NULL) {
  act <- substitute(object)
  value <- eval.parent(act)
  testthat::expect(
    is.null(value),
    sprintf("%s is not NULL", if (is.null(label)) deparse(act) else label),
    info = info
  )
  invisible(value)
}
