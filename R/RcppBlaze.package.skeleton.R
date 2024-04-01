## Copyright (C) 2010 - 2013 Dirk Eddelbuettel, Romain Francois and Douglas Bates
## Copyright (C) 2014        Dirk Eddelbuettel
## Copyright (C) 2017 - 2024 Ching-Chuan Chen
##
## This file is based on RcppArmadillo.package.skeleton.R,
## RcppArmadillo.package.skeleton.Rd, RcppEigen.package.skeleton.R
## and RcppEigen.package.skeleton.Rd from RcppArmadillo and RcppEigen.
##
## This file is part of RcppBlaze.
##
## RcppBlaze is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppBlaze is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppBlaze  If not, see <http://www.gnu.org/licenses/>.

#' Create a skeleton for a new package that intends to use RcppBlaze
#'
#' \code{RcppBlaze.package.skeleton} automates the creation of
#' a new source package that intends to use features of RcppBlaze.
#' It is based on the \link[utils]{package.skeleton} function
#' which it executes first.
#'
#' In addition to \link[utils]{package.skeleton} :
#'
#' The \samp{DESCRIPTION} file gains a Depends line requesting that
#' the package depends on Rcpp and RcppArmadillo and
#' a LinkingTo line so that the package finds Rcpp and RcppArmadillo header files.
#'
#' The \samp{NAMESPACE}, if any, gains a \code{useDynLib} directive.
#'
#' The \samp{src} directory is created if it does not exists and
#' a \samp{Makevars} file is added setting the environment variable
#' \samp{PKG_LIBS} to accomodate the necessary flags
#' to link with the Rcpp library.
#'
#' If the \code{example_code} argument is set to \code{TRUE},
#' example files \samp{rcpparma_hello_world.h} and \samp{rcpparma_hello_world.cpp}
#' are also created in the \samp{src}. An R file \samp{rcpparma_hello_world.R} is
#' expanded in the \samp{R} directory, the \code{rcpparma_hello_world} function
#' defined in this files makes use of the C++ function \samp{rcpparma_hello_world}
#' defined in the C++ file. These files are given as an example and should
#' eventually by removed from the generated package.
#'
#' @param name,list,environment,path,force,code_files see \link[utils]{package.skeleton}.
#' @param example_code If \code{TRUE}, example c++ code using \strong{RcppBlaze} is added to the package.
#' @return Nothing, used for its side effects.
#' @references
#' Read the \emph{Writing R Extensions} manual for more details.
#'
#' Once you have created a \emph{source} package you need to install it:
#' see the \emph{R Installation and Administration} manual,
#' \code{\link{INSTALL}} and \code{\link{install.packages}}.
#' @seealso \link[utils]{package.skeleton}
#' @examples
#' \dontrun{
#' RcppBlaze.package.skeleton("foobar")
#' }
#' @importFrom utils package.skeleton packageDescription
#' @export
RcppBlaze.package.skeleton <- function(
    name = "anRpackage",
    list = character(),
    environment = .GlobalEnv,
    path = ".",
    force = FALSE,
    code_files = character(),
    example_code = TRUE
) {

  env <- parent.frame(1)

  if (!length(list)) {
    fake <- TRUE
    assign("Rcpp.fake.fun", function() {}, envir = env)
  } else {
    fake <- FALSE
  }

  haveKitten <- requireNamespace("pkgKitten", quietly=TRUE)
  skelFunUsed <- ifelse(haveKitten, pkgKitten::kitten, package.skeleton)
  skelFunName <- ifelse(haveKitten, "kitten", "package.skeleton")
  message("\nCalling ", skelFunName, " to create basic package.")

  ## first let the traditional version do its business
  call <- match.call()
  call[[1]] <- skelFunUsed
  if (!haveKitten) {                 # in the package.skeleton() case
    if ("example_code" %in% names(call)) {
      call[["example_code"]] <- NULL    # remove the example_code argument
    }
    if (fake) {
      call[["list"]] <- "Rcpp.fake.fun"
    }
  }

  tryCatch(
    eval(call, envir = env),
    error = function(e) {
      cat(paste(e, "\n")) # print error
      stop(paste("error while calling `", skelFunName, "`", sep=""))
    }
  )

  message("\nAdding RcppBlaze settings")

  ## now pick things up
  root <- file.path(path, name)

  ## Add Rcpp to the DESCRIPTION
  DESCRIPTION <- file.path(root, "DESCRIPTION")
  if (file.exists(DESCRIPTION)) {
    x <- cbind(
      read.dcf(DESCRIPTION),
      "Imports" = sprintf(
        "Rcpp (>= %s), RcppBlaze (>= %s)",
        packageDescription("Rcpp")[["Version"]],
        packageDescription("BH")[["Version"]],
        packageDescription("RcppBlaze")[["Version"]]
      ),
      "LinkingTo" = "Rcpp, RcppBlaze"
    )
    write.dcf(x, file = DESCRIPTION)
    message(" >> added Imports: Rcpp, RcppBlaze")
    message(" >> added LinkingTo: Rcpp, RcppBlaze")
  }

  ## add a useDynLib to NAMESPACE,
  NAMESPACE <- file.path(root, "NAMESPACE")
  lines <- readLines(NAMESPACE)
  if (!grepl("useDynLib", lines)) {
    lines <- c(
      sprintf("useDynLib(%s)", name),
      "import(RcppBlaze)",
      "importFrom(Rcpp, evalCpp)",        ## ensures Rcpp instantiation
      lines
    )
    writeLines(lines, con = NAMESPACE)
    message(" >> added useDynLib and importFrom directives to NAMESPACE")
  }

  ## lay things out in the src directory
  src <- file.path(root, "src")
  if (!file.exists(src)) {
    dir.create(src)
  }
  man <- file.path(root, "man")
  if (!file.exists(man)) {
    dir.create(man)
  }
  skeleton <- system.file("skeleton", package = "RcppBlaze")
  Makevars <- file.path(src, "Makevars")
  if (!file.exists(Makevars)) {
    file.copy(file.path(skeleton, "Makevars"), Makevars)
    message(" >> added Makevars file with RcppBlaze settings")
  }

  if (example_code) {
    file.copy(file.path(skeleton, "rcppblaze_hello_world.cpp"), src)
    message(" >> added example src file using Blaze classes")
    file.copy(file.path(skeleton, "rcppblaze_hello_world.Rd"), man)
    message(" >> added example Rd file for using Blaze classes")

    Rcpp::compileAttributes(root)
    message(" >> invoked Rcpp::compileAttributes to create wrappers")
  }

  if (fake) {
    rm("Rcpp.fake.fun", envir = env)
    unlink(file.path(root, "R"  , "Rcpp.fake.fun.R"))
    unlink(file.path(root, "man", "Rcpp.fake.fun.Rd"))
  }

  invisible(NULL)
}
