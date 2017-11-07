#include "defs.h"
#ifndef SFC
#define SFC

#define block_size 1024

texture<float, 2, cudaReadModeElementType> texI;
texture<float, 2, cudaReadModeElementType> texPC;

#  define HM_SAFE_CALL(call) { \
    cudaError err = call;  \
    if(cudaSuccess != err) {  \
        sprintf(aD->message, "Cuda error in file '%s' in line %i: %s", __FILE__, __LINE__, cudaGetErrorString( err) ); \
		printError(aD); \
        return false;   \
    } }

#  define HMFFT_SAFE_CALL( call) {                                           \
    char * estring;                                                          \
    cufftResult err = call;                                                  \
    if( CUFFT_SUCCESS != err) {                                              \
		 if (err == CUFFT_SUCCESS) estring=strdup("CUFFT_SUCCESS");    \
		 if (err == CUFFT_INVALID_PLAN) estring=strdup("CUFFT_INVALID_PLAN");    \
		 if (err == CUFFT_ALLOC_FAILED) estring=strdup("CUFFT_ALLOC_FAILED");    \
		 if (err == CUFFT_INVALID_TYPE) estring=strdup("CUFFT_INVALID_TYPE");    \
		 if (err == CUFFT_INVALID_VALUE) estring=strdup("CUFFT_INVALID_VALUE");    \
		 if (err == CUFFT_INTERNAL_ERROR) estring=strdup("CUFFT_INTERNAL_ERROR");    \
		 if (err == CUFFT_EXEC_FAILED) estring=strdup("CUFFT_EXEC_FAILED");    \
		 if (err == CUFFT_SETUP_FAILED) estring=strdup("CUFFT_SETUP_FAILED");    \
		 if (err == CUFFT_INVALID_SIZE) estring=strdup("CUFFT_INVALID_SIZE");    \
		sprintf(aD->message, "CUFFT error in file '%s' in line %i: %s (enum = %i)", __FILE__, __LINE__,estring, err); \
		printError(aD); \
        return false;   \
    } }

#  define HMCUT_SAFE_CALL( call)                                               \
    if( CUTTrue != call) {                                                   \
		sprintf(aD->message, "Cut error in file '%s' in line %i", __FILE__, __LINE__); \
		printError(aD); \
        return false;   \
    } 

#endif

