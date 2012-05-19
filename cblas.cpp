#ifdef _WIN32

#define blasint int
#define _Complex 

#include "f2c.h"
#include "cblas.h"

int sscal_(integer *n, real *sa, real *sx, integer *incx);

int dscal_(integer *n, doublereal *da, doublereal *dx, 
	integer *incx);

int ssyr_(char *uplo, integer *n, real *alpha, real *x, 
	integer *incx, real *a, integer *lda, ftnlen uplo_len);

int dsyr_(char *uplo, integer *n, doublereal *alpha, 
	doublereal *x, integer *incx, doublereal *a, integer *lda, ftnlen 
	uplo_len);

int strsm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, real *alpha, real *a, integer *lda, real *b, 
	integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, 
	ftnlen diag_len);

int dtrsm_(char *side, char *uplo, char *transa, char *diag, 
	integer *m, integer *n, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, 
	integer *ldb, ftnlen side_len, ftnlen uplo_len, ftnlen transa_len, 
	ftnlen diag_len);

int sgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, real *alpha, real *a, integer *lda, real *b, integer *
	ldb, real *beta, real *c__, integer *ldc, ftnlen transa_len, ftnlen 
	transb_len);

int dgemm_(char *transa, char *transb, integer *m, integer *
	n, integer *k, doublereal *alpha, doublereal *a, integer *lda, doublereal *b, integer *
	ldb, doublereal *beta, doublereal *c__, integer *ldc, ftnlen transa_len, ftnlen 
	transb_len);

void cblas_sscal(blasint N, float alpha, float *X, blasint incX) {
	sscal_(&N, &alpha, X, &incX);
}

void cblas_dscal(blasint N, double alpha, double *X, blasint incX) {
	dscal_(&N, &alpha, X, &incX);
}

void cblas_ssyr(enum CBLAS_ORDER order, enum CBLAS_UPLO uplo, blasint n, float alpha, 
	float *x, blasint incx, float *a, blasint lda)
{
	ssyr_(uplo == CblasLower ? "L" : "U", &n, &alpha, x, &incx, a, &lda, 1);
}

void cblas_dsyr(enum CBLAS_ORDER order, enum CBLAS_UPLO uplo, blasint n, double alpha, 
	double *x, blasint incx, double *a, blasint lda)
{
	dsyr_(uplo == CblasLower ? "L" : "U", &n, &alpha, x, &incx, a, &lda, 1);
}

void cblas_strsm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
				 enum CBLAS_DIAG Diag, blasint M, blasint N, float alpha, float *A, blasint lda, float *B, blasint ldb) {

	 strsm_(Side == CblasLeft ? "L" : "R", Uplo == CblasLower ? "L" : "U", TransA == CblasNoTrans ? "N" : "T",
		 Diag == CblasNonUnit ? "N" : "U", &M, &N, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);	

}

void cblas_dtrsm(enum CBLAS_ORDER Order, enum CBLAS_SIDE Side, enum CBLAS_UPLO Uplo, enum CBLAS_TRANSPOSE TransA,
				 enum CBLAS_DIAG Diag, blasint M, blasint N, double alpha, double *A, blasint lda, double *B, blasint ldb) {

	 dtrsm_(Side == CblasLeft ? "L" : "R", Uplo == CblasLower ? "L" : "U", TransA == CblasNoTrans ? "N" : "T",
		 Diag == CblasNonUnit ? "N" : "U", &M, &N, &alpha, A, &lda, B, &ldb, 1, 1, 1, 1);	

}

void cblas_sgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint M, blasint N, blasint K,
				 float alpha, float *A, blasint lda, float *B, blasint ldb, float beta, float *C, blasint ldc) {

	sgemm_(TransA == CblasNoTrans ? "N" : "T", TransB == CblasNoTrans ? "N" : "T", 
		&M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);

}

void cblas_dgemm(enum CBLAS_ORDER Order, enum CBLAS_TRANSPOSE TransA, enum CBLAS_TRANSPOSE TransB, blasint M, blasint N, blasint K,
				 double alpha, double *A, blasint lda, double *B, blasint ldb, double beta, double *C, blasint ldc) {

	dgemm_(TransA == CblasNoTrans ? "N" : "T", TransB == CblasNoTrans ? "N" : "T", 
		&M, &N, &K, &alpha, A, &lda, B, &ldb, &beta, C, &ldc, 1, 1);

}

#endif