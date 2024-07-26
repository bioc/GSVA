#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* prototypes of functions to be registered */

SEXP
matrix_density_R(double* X, double* Y, int* n_density_samples,
                 int* n_test_samples, int* n_genes, int* rnaseq);

SEXP
ks_matrix_R(SEXP XR, SEXP sidxsR, SEXP n_genesR, SEXP geneset_idxsR,
            SEXP n_genesetR, SEXP tauR, SEXP n_samplesR, SEXP mx_diffR, SEXP abs_rnkR);

SEXP
ecdfvals_sparse_to_sparse_R(SEXP XCspR, SEXP XRspR);

SEXP
ecdfvals_sparse_to_dense_R(SEXP XCspR, SEXP XRspR);

SEXP
ecdfvals_dense_to_dense_R(SEXP XR);

/* registration of C-entry points */

static R_CallMethodDef callMethods[] = {
  {"ks_matrix_R", (DL_FUNC) &ks_matrix_R, 9},
  {"matrix_density_R", (DL_FUNC) &matrix_density_R, 6},
  {"ecdfvals_sparse_to_sparse_R", (DL_FUNC) &ecdfvals_sparse_to_sparse_R, 2},
  {"ecdfvals_sparse_to_dense_R", (DL_FUNC) &ecdfvals_sparse_to_dense_R, 2},
  {"ecdfvals_dense_to_dense_R", (DL_FUNC) &ecdfvals_dense_to_dense_R, 1},
  {NULL, NULL, 0}
};

/* global variables */
SEXP Matrix_DimNamesSym,
     Matrix_DimSym,
     Matrix_xSym,
     Matrix_iSym,
     Matrix_jSym,
     Matrix_pSym;

void
R_init_GSVA(DllInfo *info) {

  R_registerRoutines(info, NULL, callMethods, NULL, NULL);

  /* from the Matrix package init.c */
  Matrix_DimNamesSym = install("Dimnames");
  Matrix_DimSym = install("Dim");
  Matrix_xSym = install("x");
  Matrix_iSym = install("i");
  Matrix_jSym = install("j");
  Matrix_pSym = install("p");

  R_useDynamicSymbols(info, TRUE);

}
