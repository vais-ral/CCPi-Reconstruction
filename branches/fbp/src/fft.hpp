
#ifndef FFT_WRAPPERS
#define FFT_WRAPPERS

template <class real_type> void fft_1d_forward(std::complex<real_type> data[],
					       const sl_int n, const int ntr)
{
  std::cerr << "Forward 1D FFT not available\n";
}

template <class real_type> void fft_1d_inverse(std::complex<real_type> data[],
					       const sl_int n, const int ntr)
{
  std::cerr << "Inverse 1D FFT not available\n";
}

#  if defined(MKL_ILP64)
#    include <mkl_dfti.h>

template <> inline void fft_1d_forward(std::complex<float> data[],
				       const sl_int n, const int ntr)
{
  DFTI_DESCRIPTOR_HANDLE desc;
  (void) DftiCreateDescriptor(&desc, DFTI_SINGLE, DFTI_COMPLEX, 1, n);
  DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, ntr);
  DftiSetValue(desc, DFTI_INPUT_DISTANCE, n);
  DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, n);
  DftiCommitDescriptor(desc);
  DftiComputeForward(desc, data);
  DftiFreeDescriptor(&desc);
}

template <> inline void fft_1d_forward(std::complex<double> data[],
				       const sl_int n, const int ntr)
{
  DFTI_DESCRIPTOR_HANDLE desc;
  (void) DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, n);
  DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, ntr);
  DftiSetValue(desc, DFTI_INPUT_DISTANCE, n);
  DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, n);
  DftiCommitDescriptor(desc);
  DftiComputeForward(desc, data);
  DftiFreeDescriptor(&desc);
}

template <> inline void fft_1d_inverse(std::complex<float> data[],
				       const sl_int n, const int ntr)
{
  DFTI_DESCRIPTOR_HANDLE desc;
  (void) DftiCreateDescriptor(&desc, DFTI_SINGLE, DFTI_COMPLEX, 1, n);
  DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, ntr);
  DftiSetValue(desc, DFTI_INPUT_DISTANCE, n);
  DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, n);
  DftiCommitDescriptor(desc);
  DftiComputeBackward(desc, data);
  DftiFreeDescriptor(&desc);
}

template <> inline void fft_1d_inverse(std::complex<double> data[],
				       const sl_int n, const int ntr)
{
  DFTI_DESCRIPTOR_HANDLE desc;
  (void) DftiCreateDescriptor(&desc, DFTI_DOUBLE, DFTI_COMPLEX, 1, n);
  DftiSetValue(desc, DFTI_NUMBER_OF_TRANSFORMS, ntr);
  DftiSetValue(desc, DFTI_INPUT_DISTANCE, n);
  DftiSetValue(desc, DFTI_OUTPUT_DISTANCE, n);
  DftiCommitDescriptor(desc);
  DftiComputeBackward(desc, data);
  DftiFreeDescriptor(&desc);
}

#  elif defined(USE_FFTW)
#    error "FFTW - todo"
#  endif // FFTW

#endif // FFT_WRAPPERS
