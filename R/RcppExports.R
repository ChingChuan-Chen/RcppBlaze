# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' The version of Blaze used in RcppBlaze
#'
#' To return the version of Blaze used in RcppBlaze.
#'
#' @param single A logical value indicates which type to return. If TRUE, it returns an integer. If FALSE, it returns a named vector.
#' @return A number or a named vector to represent the version of \code{blaze} depending on the input, \code{single}.
#' @seealso Blaze header file \code{blaze/system/Version.h}.
#' @examples
#' blaze_version(FALSE)
#' @export
blaze_version <- function(single) {
    .Call(`_RcppBlaze_blaze_version`, single)
}

#' Set/Get the random number generator for blaze with given seed
#'
#' @param seed A positive integer to specify the seed value for the random number generator.
#' @return No return value.
#' @rdname blaze_seed
#' @export
blaze_set_seed <- function(seed) {
    invisible(.Call(`_RcppBlaze_blaze_set_seed`, seed))
}

#' @rdname blaze_seed
#' @export
blaze_get_seed <- function() {
    .Call(`_RcppBlaze_blaze_get_seed`)
}

#' Set/Get the Number of Threads used in blaze
#'
#' @param n The number of threads to set in blaze.
#' @return \code{blaze_get_threads} returns an integer and \code{blaze_set_threads} returns nothing.
#' @seealso blaze wiki: \url{https://bitbucket.org/blaze-lib/blaze/wiki/Shared\%20Memory\%20Parallelization}.
#' @rdname blaze_threads
#' @export
blaze_set_num_threads <- function(n) {
    invisible(.Call(`_RcppBlaze_blaze_set_num_threads`, n))
}

#' @rdname blaze_threads
#' @export
blaze_get_num_threads <- function() {
    .Call(`_RcppBlaze_blaze_get_num_threads`)
}

#' linear model fitting function based on RcppBlaze
#'
#' \code{fastLmPure} provides the estimates of the linear model based on \strong{RcppBlaze}.
#'
#' \code{fastLm} estimates the linear model using the \code{solve}.
#'
#' @param X A model matrix.
#' @param y A response vector.
#' @param type A integer. 0 is QR solver, 1 is LDLT solver, 2 is LLT sovler and 3 is LU solver.
#' @return A list containing coefficients, standard errors, rank of model matrix,
#'   degree of freedom of residuals, residuals, the standard deviation of random errors and
#'   fitted values.
#' @examples
#' # according to fastLm example in RcppArmadillo
#' data(trees, package="datasets")
#' flm <- fastLmPure(cbind(1, log(trees$Girth)), log(trees$Volume), 0)
#' print(flm)
#' @export
fastLmPure <- function(X, y, type) {
    .Call(`_RcppBlaze_fastLmPure`, X, y, type)
}

