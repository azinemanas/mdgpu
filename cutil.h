#ifndef _CUTIL_H
#define _CUTIL_H

#include "types.h"

#define cutilSafeCall(err)           __cudaSafeCall      (err, __FILE__, __LINE__)
#define cutilCheckMsg(msg)           __cutilCheckMsg     (msg, __FILE__, __LINE__)

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
    cudaThreadSynchronize();
    if( cudaSuccess != err) {
  printf("%s(%i) : cudaSafeCall() Runtime API error (%i): %s.\n",
                file, line, err, cudaGetErrorString( err) );
        exit(-1);
    }
}

inline void __cutilCheckMsg( const char *errorMessage, const char *file, const int line )
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        printf("%s(%i) : cutilCheckMsg() CUTIL CUDA error : %s : %s.\n",
                file, line, errorMessage, cudaGetErrorString( err) );
        exit(-1);
    }
//#ifdef _DEBUG
    err = cudaThreadSynchronize();
    if( cudaSuccess != err) {
  printf("%s(%i) : cutilCheckMsg cudaThreadSynchronize error: %s : %s.\n",
                file, line, errorMessage, cudaGetErrorString( err) );
        exit(-1);
    }
//#endif
}

#endif