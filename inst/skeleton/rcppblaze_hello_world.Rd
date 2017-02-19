\name{RcppBlaze-Functions}
\alias{rcppblaze_hello_world}
\alias{rcppblaze_innerproduct}
\alias{rcppblaze_outerproduct}
\alias{rcppblaze_bothproducts}
\title{Set of functions in example RcppBlaze package}
\description{
  These four functions are created when
  \code{RcppBlaze.package.skeleton()} is invoked to create a
  skeleton packages.
}
\usage{
rcppblaze_hello_world()
rcppblaze_outerproduct(x)
rcppblaze_innerproduct(x)
rcppblaze_bothproducts(x)
}
\arguments{
  \item{x}{a numeric vector}
}
\value{
  \code{rcppblaze_hello_world()} does not return a value, but displays a
  message to the console.

  \code{rcppblaze_outerproduct()} returns a numeric matrix computed as the
  outer (vector) product of \code{x}.

  \code{rcppblaze_innerproduct()} returns a double computer as the inner
  (vector) product of \code{x}.

  \code{rcppblaze_bothproducts()} returns a list with both the outer and
  inner products.

}
\details{
  These are example functions which should be largely
  self-explanatory. Their main benefit is to demonstrate how to write a
  function using the Blaze C++ classes, and to have to such a
  function accessible from R.
}
\references{
  See the documentation for Blaze, and RcppBlaze, for more details.
}
\examples{
  x <- sqrt(1:4)
  rcppblaze_innerproduct(x)
  rcppblaze_outerproduct(x)
}
\author{Dirk Eddelbuettel}

