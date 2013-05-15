
#ifndef BLAS_WRAPPERS
#define BLAS_WRAPPERS

#ifndef MKL

/* Only implements the simplest version of ddot which is used in this code 
   here, with incrx and incry = 1 hardcoded */
template <class real_type>
real_type ddot(const long n, const real_type x[], const int incx,
	       const real_type y[], const int incy)
{
  real_type sum = 0.0;
  for (long i = 0; i < n; i++)
    sum += x[i] * y[i];
  return sum;
}

template <class real_type>
real_type dnrm2(const long n, const real_type x[], const int incx)
{
  real_type nrm2 = 0.0;
  for (long i = 0; i < n; i++)
    nrm2 += x[i] * x[i];
  return sqrt(nrm2);
}

/* Only implements the simplest version of daxpy which is used in this code 
   here, with incrx and incry = 1 hardcoded. *y points at the results
   (overwrite) */
template <class real_type>
void daxpy(const long n, const real alpha, const real_type x[],
	   const int incx, real_type y[], const int incy)
{
  real_type a = real_type(alpha);
  for (long i = 0; i < n; i++)
    y[i] += a * x[i];
}

template <class real_type>
void dcopy(const long n, const real_type x[], const int incx, real_type y[],
	   const int incy)
{
  for (long i = 0; i < n; i++)
    y[i] = x[i];
}

#endif // MKL

#endif // BLAS_WRAPPERS
