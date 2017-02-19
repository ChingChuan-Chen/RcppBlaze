// [[Rcpp::depends(RcppBlaze)]]
#include <RcppBlaze.h>

// [[Rcpp::export]]
blaze::DenseColumn< blaze::DynamicMatrix<double> > test_column_dbl(blaze::DynamicMatrix<double> x, int idx) {
  if ( idx < 0 || idx >= (int) x.columns())
    Rcpp::stop("idx does not fit the dimension of x.");
  blaze::DenseColumn< blaze::DynamicMatrix<double> > y = blaze::column( x, (size_t) idx );
  return y;
}

// [[Rcpp::export]]
blaze::DenseColumn< blaze::DynamicMatrix< std::complex<double> > > test_column_cpl(blaze::DynamicMatrix< std::complex<double> > x, int idx) {
  if ( idx < 0 || idx >= (int) x.columns())
    Rcpp::stop("idx does not fit the dimension of x.");
  blaze::DenseColumn< blaze::DynamicMatrix< std::complex<double> > > y = blaze::column( x, (size_t) idx );
  return y;
}

// [[Rcpp::export]]
blaze::DenseRow< blaze::DynamicMatrix<double> > test_row_dbl(blaze::DynamicMatrix<double> x, int idx) {
  if ( idx < 0 || idx >= (int) x.rows())
    Rcpp::stop("idx does not fit the dimension of x.");
  blaze::DenseRow< blaze::DynamicMatrix<double> > y = blaze::row( x, (size_t) idx );
  return y;
}

// [[Rcpp::export]]
blaze::DenseRow< blaze::DynamicMatrix< std::complex<double> > > test_row_cpl(blaze::DynamicMatrix< std::complex<double> > x, int idx) {
  if ( idx < 0 || idx >= (int) x.rows())
    Rcpp::stop("idx does not fit the dimension of x.");
  blaze::DenseRow< blaze::DynamicMatrix< std::complex<double> > > y = blaze::row( x, (size_t) idx );
  return y;
}

// [[Rcpp::export]]
blaze::DenseSubvector< blaze::DynamicVector<double> > test_subvec_col_dbl(blaze::DynamicVector<double> x, int x1, int x2 ) {
  if ( x1 < 0 || x2 >= (int) x.size() )
    Rcpp::stop("x1 and x2 does not fit the dimension of x.");
  blaze::DenseSubvector< blaze::DynamicVector<double> > y = blaze::subvector( x, (size_t) x1, (size_t) x2 );
  return y;
}

// [[Rcpp::export]]
blaze::DenseSubvector< blaze::DynamicVector<double, blaze::rowVector> > test_subvec_row_dbl(blaze::DynamicVector<double, blaze::rowVector> x, int x1, int x2 ) {
  if ( x1 < 0 || x2 >= (int) x.size() )
    Rcpp::stop("x1 and x2 does not fit the dimension of x.");
  blaze::DenseSubvector< blaze::DynamicVector<double, blaze::rowVector> > y = blaze::subvector( x, (size_t) x1, (size_t) x2 );
  return y;
}

// [[Rcpp::export]]
blaze::DenseSubvector< blaze::DynamicVector< std::complex<double> > > test_subvec_col_cpl(blaze::DynamicVector< std::complex<double> > x, int x1, int x2) {
  if ( x1 < 0 || x2 >= (int) x.size() )
    Rcpp::stop("x1 and x2 does not fit the dimension of x.");
  blaze::DenseSubvector< blaze::DynamicVector< std::complex<double> > > y = blaze::subvector( x, (size_t) x1, (size_t) x2 );
  return y;
}

// [[Rcpp::export]]
blaze::DenseSubvector< blaze::DynamicVector< std::complex<double>, blaze::rowVector > > test_subvec_row_cpl(blaze::DynamicVector< std::complex<double>, blaze::rowVector > x, int x1, int x2) {
  if ( x1 < 0 || x2 >= (int) x.size() )
    Rcpp::stop("x1 and x2 does not fit the dimension of x.");
  blaze::DenseSubvector< blaze::DynamicVector< std::complex<double>, blaze::rowVector > > y = blaze::subvector( x, (size_t) x1, (size_t) x2 );
  return y;
}

// [[Rcpp::export]]
blaze::DenseSubmatrix< blaze::DynamicMatrix<double> > test_submat_dbl(blaze::DynamicMatrix<double> x, int x1, int x2, int y1, int y2 ) {
  if ( x1 < 0 || x2 >= (int) x.rows() || y1 < 0 || y2 >= (int) x.columns() )
    Rcpp::stop("x1, x2, y1 and y2 does not fit the dimension of x.");
  blaze::DenseSubmatrix< blaze::DynamicMatrix<double> > y = blaze::submatrix( x, (size_t) x1, (size_t) x2, (size_t) y1, (size_t) y2 );
  return y;
}

// [[Rcpp::export]]
blaze::DenseSubmatrix< blaze::DynamicMatrix< std::complex<double> > > test_submat_cpl(blaze::DynamicMatrix< std::complex<double> > x, int x1, int x2, int y1, int y2 ) {
  if ( x1 < 0 || x2 >= (int) x.rows() ||  y1 < 0 || y2 >= (int) x.columns() )
    Rcpp::stop("x1, x2, y1 and y2 does not fit the dimension of x.");
  blaze::DenseSubmatrix< blaze::DynamicMatrix< std::complex<double> > > y = blaze::submatrix( x, (size_t) x1, (size_t) x2, (size_t) y1, (size_t) y2 );
  return y;
}

// [[Rcpp::export]]
blaze::SparseColumn< blaze::CompressedMatrix<double> > test_sp_column_dbl(blaze::CompressedMatrix<double> x, int idx) {
  if ( idx < 0 || idx >= (int) x.columns())
    Rcpp::stop("idx does not fit the dimension of x.");
  blaze::SparseColumn< blaze::CompressedMatrix<double> > y = blaze::column( x, (size_t) idx );
  return y;
}

// [[Rcpp::export]]
blaze::SparseRow< blaze::CompressedMatrix<double> > test_sp_row_dbl(blaze::CompressedMatrix<double> x, int idx) {
  if ( idx < 0 || idx >= (int) x.rows())
    Rcpp::stop("idx does not fit the dimension of x.");
  blaze::SparseRow< blaze::CompressedMatrix<double> > y = blaze::row( x, (size_t) idx );
  return y;
}

// [[Rcpp::export]]
blaze::SparseSubvector< blaze::CompressedVector<double> > test_sp_subcolvec_dbl(blaze::CompressedVector<double> x, int x1, int x2) {
  if ( x1 < 0 || x2 >= (int) x.size())
    Rcpp::stop("x1 and x2 does not fit the dimension of x.");
  blaze::SparseSubvector< blaze::CompressedVector<double> > y = blaze::subvector( x, (size_t) x1, (size_t) x2 );
  return y;
}

// [[Rcpp::export]]
blaze::SparseSubvector< blaze::CompressedVector<double,blaze::rowVector> > test_sp_subrowvec_dbl(blaze::CompressedVector<double,blaze::rowVector> x, int x1, int x2) {
  if ( x1 < 0 || x2 >= (int) x.size())
    Rcpp::stop("x1 and x2 does not fit the dimension of x.");
  blaze::SparseSubvector< blaze::CompressedVector<double,blaze::rowVector> > y = blaze::subvector( x, (size_t) x1, (size_t) x2 );
  return y;
}

// [[Rcpp::export]]
blaze::SparseSubmatrix< blaze::CompressedMatrix<double> > test_sp_submat_dbl(blaze::CompressedMatrix<double> x, int x1, int x2, int y1, int y2 ) {
  if ( x1 < 0 || x2 >= (int) x.rows() || y1 < 0 || y2 >= (int) x.columns() )
    Rcpp::stop("x1, x2, y1 and y2 does not fit the dimension of x.");
  blaze::SparseSubmatrix< blaze::CompressedMatrix<double> > y = blaze::submatrix( x, (size_t) x1, (size_t) x2, (size_t) y1, (size_t) y2 );
  return y;
}
