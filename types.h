#ifndef _TYPES_H
#define _TYPES_H

#if 1
	#define FLOTANTE float
	#define mdgpu_cublasXsyr cublasSsyr
	#define mdgpu_cublasXtrsm cublasStrsm
	#define mdgpu_cublasXgemm cublasSgemm

	#define mdgpu_cblasXscal cblas_sscal
	#define mdgpu_cblasXsyr cblas_ssyr
	#define mdgpu_cblasXtrsm cblas_strsm
	#define mdgpu_cblasXgemm cblas_sgemm
#else
	#define FLOTANTE double
	#define mdgpu_cublasXsyr cublasDsyr
	#define mdgpu_cublasXtrsm cublasDtrsm
	#define mdgpu_cublasXgemm cublasDgemm

	#define mdgpu_cblasXscal cblas_dscal
	#define mdgpu_cblasXsyr cblas_dsyr
	#define mdgpu_cblasXtrsm cblas_dtrsm
	#define mdgpu_cblasXgemm cblas_dgemm
#endif

#endif
