## Copyright (C) 2017 Chingchuan Chen
## Earlier copyrights Dirk Eddelbuettel, Romain Francois and Douglas Bates,
##   Gregor Gorjanc, Martin Maechler and Murray Stokely as detailed below
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
## along with RcppBlaze.  If not, see <http://www.gnu.org/licenses/>.

## doRUnit.R --- Run RUnit tests
##
## with credits to package fUtilities in RMetrics
## which credits Gregor Gojanc's example in CRAN package  'gdata'
## as per the R Wiki http://wiki.r-project.org/rwiki/doku.php?id=developers:runit
## and changed further by Martin Maechler
## and more changes by Murray Stokely in HistogramTools
## and then used adapted in RProtoBuf
## and now used in Rcpp, RcppArmadillo and Here
##
## Chingchuan Chen, Feb 2017

if (require("RUnit", quietly = TRUE)) {
  pkg <- "RcppBlaze"
  require(pkg, character.only = TRUE)
  path <- system.file("unitTests", package = pkg)
  stopifnot(file.exists(path), file.info(path.expand(path))$isdir)

  ## without this, we get unit test failures
  Sys.setenv(R_TESTS = "")
  source(file.path(path, "runTests.R"), echo = TRUE)
} else {
  print( "package RUnit not available, cannot run unit tests" )
}
