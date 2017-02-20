//[[Rcpp::depends(RcppBlaze)]]
#include <RcppBlaze.h>

// [[Rcpp::export]]
blaze::DiagonalMatrix< blaze::DynamicMatrix<double> > test_diag_matrix_dbl(blaze::DynamicMatrix<double> x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::DiagonalMatrix< blaze::DynamicMatrix<double> > A( x.rows(), x.columns() );
  for ( size_t i=0; i<x.rows(); ++i)
    A(i, i) = x(i, i);
  return A;
}

// [[Rcpp::export]]
blaze::DiagonalMatrix< blaze::DynamicMatrix< std::complex<double> > > test_diag_matrix_cpl(blaze::DynamicMatrix< std::complex<double> > x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::DiagonalMatrix< blaze::DynamicMatrix< std::complex<double> > > A( x.rows(), x.columns() );
  for ( size_t i=0; i<x.rows(); ++i)
    A(i, i) = x(i, i);
  return A;
}

// [[Rcpp::export]]
blaze::LowerMatrix< blaze::DynamicMatrix<double> > test_lower_matrix_dbl(blaze::DynamicMatrix<double> x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::LowerMatrix< blaze::DynamicMatrix<double> > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=j; i<x.rows(); ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::LowerMatrix< blaze::DynamicMatrix< std::complex<double> > > test_lower_matrix_cpl(blaze::DynamicMatrix< std::complex<double> > x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::LowerMatrix< blaze::DynamicMatrix< std::complex<double> > > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=j; i<x.rows(); ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::UpperMatrix< blaze::DynamicMatrix<double> > test_upper_matrix_dbl(blaze::DynamicMatrix<double> x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::UpperMatrix< blaze::DynamicMatrix<double> > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=0; i<=j; ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::UpperMatrix< blaze::DynamicMatrix< std::complex<double> > > test_upper_matrix_cpl(blaze::DynamicMatrix< std::complex<double> > x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::UpperMatrix< blaze::DynamicMatrix< std::complex<double> > > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=0; i<=j; ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::StrictlyLowerMatrix< blaze::DynamicMatrix<double> > test_strictlylower_matrix_dbl(blaze::DynamicMatrix<double> x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::StrictlyLowerMatrix< blaze::DynamicMatrix<double> > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=j+1; i<x.rows(); ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::StrictlyLowerMatrix< blaze::DynamicMatrix< std::complex<double> > > test_strictlylower_matrix_cpl(blaze::DynamicMatrix< std::complex<double> > x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::StrictlyLowerMatrix< blaze::DynamicMatrix< std::complex<double> > > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=j+1; i<x.rows(); ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::StrictlyUpperMatrix< blaze::DynamicMatrix<double> > test_strictlyupper_matrix_dbl(blaze::DynamicMatrix<double> x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::StrictlyUpperMatrix< blaze::DynamicMatrix<double> > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=0; i<j; ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::StrictlyUpperMatrix< blaze::DynamicMatrix< std::complex<double> > > test_strictlyupper_matrix_cpl(blaze::DynamicMatrix< std::complex<double> > x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::StrictlyUpperMatrix< blaze::DynamicMatrix< std::complex<double> > > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=0; i<j; ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::UniLowerMatrix< blaze::DynamicMatrix<double> > test_unilower_matrix_dbl(blaze::DynamicMatrix<double> x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::UniLowerMatrix< blaze::DynamicMatrix<double> > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=j+1; i<x.rows(); ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::UniLowerMatrix< blaze::DynamicMatrix< std::complex<double> > > test_unilower_matrix_cpl(blaze::DynamicMatrix< std::complex<double> > x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::UniLowerMatrix< blaze::DynamicMatrix< std::complex<double> > > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=j+1; i<x.rows(); ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::UniUpperMatrix< blaze::DynamicMatrix<double> > test_uniupper_matrix_dbl(blaze::DynamicMatrix<double> x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::UniUpperMatrix< blaze::DynamicMatrix<double> > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=0; i<j; ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::UniUpperMatrix< blaze::DynamicMatrix< std::complex<double> > > test_uniupper_matrix_cpl(blaze::DynamicMatrix< std::complex<double> > x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::UniUpperMatrix< blaze::DynamicMatrix< std::complex<double> > > A( x.rows(), x.columns() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=0; i<j; ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::SymmetricMatrix< blaze::DynamicMatrix<double> > test_symmetric_matrix( blaze::DynamicMatrix<double> x) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::SymmetricMatrix< blaze::DynamicMatrix<double> > A( x.rows() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=0; i<=j; ++i)
      A(i, j) = x(i, j);
  return A;
}

// [[Rcpp::export]]
blaze::HermitianMatrix< blaze::DynamicMatrix< std::complex<double> > > test_hermitian_matrix( blaze::DynamicMatrix< std::complex<double> > x ) {
  if ( x.rows() != x.columns() )
    Rcpp::stop("x must be square matrix.");
  blaze::HermitianMatrix< blaze::DynamicMatrix< std::complex<double> > > A( x.rows() );
  for ( size_t j=0; j<x.columns(); ++j)
    for ( size_t i=0; i<=j; ++i)
      A(i, j) = x(i, j);
  return A;
}

