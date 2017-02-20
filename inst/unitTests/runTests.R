## Copyright (C)       2017 Chingchuan Chen
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

# This script is modified from 'inst/unitTests/runTests.R' in Rcpp, RcppArmadillo and RcppEigen.
require(Rcpp)
require(Matrix)
pkg <- "RcppBlaze"

if (require("RUnit", quietly = TRUE)) {

  library(package = pkg, character.only = TRUE)
  if (!(exists("path") && file.exists(path)))
    path <- system.file("unitTests", package = pkg)

  ## --- Testing ---

  ## Define tests
  testSuite <- defineTestSuite(paste(pkg, "unit testing"), path, "^rUnitTest-.+\\.[rR]$")

  if (interactive()) {
    cat("Now have RUnit Test Suite 'testSuite' for package '", pkg, "' :\n", sep='')
    str(testSuite)
    cat('', "Consider doing",
        "\t  tests <- runTestSuite(testSuite)", "\nand later",
        "\t  printTextProtocol(tests)", '', sep="\n")
  } else { ## run from shell / Rscript / R CMD Batch / ...
    ## Run
    tests <- runTestSuite(testSuite)

    tmpDir <- tempdir()
    pathReport <- file.path(tmpDir, paste0(pkg, "_unitTests_report"))
    cat("RUnit reports are written into ", tmpDir, "/", pkg, "_unitTests_report(txt|html)", sep = "")

    ## Print results
    printTextProtocol(tests, paste0(pathReport, ".txt"))
    ## Print HTML version to a file
    ## printHTMLProtocol has problems on Mac OS X
    if (Sys.info()["sysname"] != "Darwin")
      printHTMLProtocol(tests, paste0(pathReport, ".html"))

    ## stop() if there are any failures i.e. FALSE to unit test.
    ## This will cause R CMD check to return error and stop
    if (getErrors(tests)$nFail > 0)
      stop("one of the unit tests failed")
  }
} else {
  cat("R package 'RUnit' cannot be loaded -- no unit tests run\n", "for package", pkg, "\n")
}
