template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatAddExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatKronExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, typename OP, bool SO> SEXP wrap(const blaze::DMatDMatMapExpr<MT1, MT2, OP, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::DMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSchurExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSolveExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatDMatSubExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT> SEXP wrap(const blaze::DMatDVecMultExpr<MT, VT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT, bool TF> SEXP wrap(const blaze::DMatDVecSolveExpr<MT, VT, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename OP, bool SO> SEXP wrap(const blaze::DMatMapExpr<MT, OP, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::DMatScalarDivExpr<MT, ST, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::DMatScalarMultExpr<MT, ST, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatAddExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatKronExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::DMatSMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::DMatSMatSchurExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::DMatSMatSubExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT> SEXP wrap(const blaze::DMatSVecMultExpr<MT, VT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecAddExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecCrossExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecDivExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecKronExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, typename OP, bool TF> SEXP wrap(const blaze::DVecDVecMapExpr<VT1, VT2, OP, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecMultExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, typename OP> SEXP wrap(const blaze::DVecDVecOuterExpr<VT1, VT2, OP>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecDVecSubExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename OP, bool TF> SEXP wrap(const blaze::DVecMapExpr<VT, OP, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::DVecScalarDivExpr<VT, ST, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::DVecScalarMultExpr<VT, ST, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecAddExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecCrossExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecKronExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecMultExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2> SEXP wrap(const blaze::DVecSVecOuterExpr<VT1, VT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::DVecSVecSubExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::SMatDMatKronExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::SMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatDMatSchurExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SO> SEXP wrap(const blaze::SMatDMatSubExpr<MT1, MT2, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT> SEXP wrap(const blaze::SMatDVecMultExpr<MT, VT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename OP, bool SO> SEXP wrap(const blaze::SMatMapExpr<MT, OP, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::SMatScalarDivExpr<MT, ST, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename ST, bool SO> SEXP wrap(const blaze::SMatScalarMultExpr<MT, ST, SO>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatAddExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatKronExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatMultExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatSchurExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::SMatSMatSubExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT> SEXP wrap(const blaze::SMatSVecMultExpr<MT, VT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecCrossExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecDivExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecKronExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecMultExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2> SEXP wrap(const blaze::SVecDVecOuterExpr<VT1, VT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecDVecSubExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename OP, bool TF> SEXP wrap(const blaze::SVecMapExpr<VT, OP, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::SVecScalarDivExpr<VT, ST, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename ST, bool TF> SEXP wrap(const blaze::SVecScalarMultExpr<VT, ST, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecAddExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecCrossExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecKronExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecMultExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2> SEXP wrap(const blaze::SVecSVecOuterExpr<VT1, VT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT1, typename VT2, bool TF> SEXP wrap(const blaze::SVecSVecSubExpr<VT1, VT2, TF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TDMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT> SEXP wrap(const blaze::TDMatDVecMultExpr<MT, VT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::TDMatSMatAddExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TDMatSMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::TDMatSMatSubExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT> SEXP wrap(const blaze::TDMatSVecMultExpr<MT, VT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename MT> SEXP wrap(const blaze::TDVecDMatMultExpr<VT, MT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename MT> SEXP wrap(const blaze::TDVecSMatMultExpr<VT, MT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2, bool SF, bool HF, bool LF, bool UF> SEXP wrap(const blaze::TSMatDMatMultExpr<MT1, MT2, SF, HF, LF, UF>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatDMatSchurExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatDMatSubExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT> SEXP wrap(const blaze::TSMatDVecMultExpr<MT, VT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatKronExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatMultExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatSchurExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT1, typename MT2> SEXP wrap(const blaze::TSMatSMatSubExpr<MT1, MT2>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename MT, typename VT> SEXP wrap(const blaze::TSMatSVecMultExpr<MT, VT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename MT> SEXP wrap(const blaze::TSVecDMatMultExpr<VT, MT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

template <typename VT, typename MT> SEXP wrap(const blaze::TSVecSMatMultExpr<VT, MT>& x) {
  return RcppBlaze::blaze_wrap(x);
};

