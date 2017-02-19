// [[Rcpp::depends(RcppBlaze)]]
#include <RcppBlaze.h>

// [[Rcpp::export]]
blaze::StaticMatrix< double,4UL,4UL > test_StaticMatrix_dbl_dim44( blaze::StaticMatrix< double,4UL,4UL > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::StaticMatrix< std::complex<double>,4UL,4UL > test_StaticMatrix_cpl_dim44( blaze::StaticMatrix< std::complex<double>,4UL,4UL > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::HybridMatrix< double,4UL,4UL > test_HybridMatrix_dbl_dim44( blaze::HybridMatrix< double,4UL,4UL > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::HybridMatrix< std::complex<double>,4UL,4UL > test_HybridMatrix_cpl_dim44( blaze::HybridMatrix< std::complex<double>,4UL,4UL > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::DynamicMatrix<double> test_DynamicMatrix_dbl(blaze::DynamicMatrix<double> x) {
  return x;
}

// [[Rcpp::export]]
blaze::DynamicMatrix< std::complex<double> > test_DynamicMatrix_cpl( blaze::DynamicMatrix< std::complex<double> > x ) {
  return x;
}

// [[Rcpp::export]]
blaze::CustomMatrix< double, blaze::unaligned, blaze::unpadded >
  test_CustomMatrix1_dbl( blaze::CustomMatrix< double, blaze::unaligned, blaze::unpadded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomMatrix< double, blaze::aligned, blaze::unpadded >
  test_CustomMatrix2_dbl( blaze::CustomMatrix< double, blaze::aligned, blaze::unpadded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomMatrix< double, blaze::unaligned, blaze::padded >
  test_CustomMatrix3_dbl( blaze::CustomMatrix< double, blaze::unaligned, blaze::padded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomMatrix< double, blaze::aligned, blaze::padded >
  test_CustomMatrix4_dbl( blaze::CustomMatrix< double, blaze::aligned, blaze::padded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomMatrix< std::complex<double>, blaze::unaligned, blaze::unpadded >
  test_CustomMatrix1_cpl( blaze::CustomMatrix< std::complex<double>, blaze::unaligned, blaze::unpadded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomMatrix< std::complex<double>, blaze::aligned, blaze::unpadded >
  test_CustomMatrix2_cpl( blaze::CustomMatrix< std::complex<double>, blaze::aligned, blaze::unpadded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomMatrix< std::complex<double>, blaze::unaligned, blaze::padded >
  test_CustomMatrix3_cpl( blaze::CustomMatrix< std::complex<double>, blaze::unaligned, blaze::padded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CustomMatrix< std::complex<double>, blaze::aligned, blaze::padded >
  test_CustomMatrix4_cpl( blaze::CustomMatrix< std::complex<double>, blaze::aligned, blaze::padded > x ) {
    return x;
}

// [[Rcpp::export]]
blaze::CompressedMatrix<double> test_CompressedMatrix_dbl_dgC(blaze::CompressedMatrix<double> x) {
  return x;
}

// [[Rcpp::export]]
blaze::CompressedMatrix<double, blaze::rowMajor> test_CompressedMatrix_dbl_dgR(blaze::CompressedMatrix<double, blaze::rowMajor> x) {
  return x;
}
