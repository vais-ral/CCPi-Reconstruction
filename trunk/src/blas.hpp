
#ifndef BLAS_WRAPPERS
#define BLAS_WRAPPERS

/* Only implements the simplest version of ddot which is used in this code 
   here, with incrx and incry = 1 hardcoded */
template <class real_type>
real_type ddot(const sl_int n, const real_type x[], const int incx,
	       const real_type y[], const int incy)
{
  real_type sum = 0.0;
  for (sl_int i = 0; i < n; i++)
    sum += x[i] * y[i];
  return sum;
}

template <class real_type>
real_type ddot(const sl_int n, const boost::multi_array_ref<real_type, 3> &x,
	       const int incx, const boost::multi_array_ref<real_type, 3> &y,
	       const int incy)
{
  real_type sum = 0.0;
  for (sl_int i = 0; i < n; i++)
    sum += (x.data())[i] * (y.data())[i];
  return sum;
}

template <class real_type>
real_type dnrm2(const sl_int n, const real_type x[], const int incx)
{
  real_type nrm2 = 0.0;
  for (sl_int i = 0; i < n; i++)
    nrm2 += x[i] * x[i];
  return sqrt(nrm2);
}

template <class real_type>
real_type dnrm2(const sl_int n, const boost::multi_array_ref<real_type, 3> &x,
		const int incx)
{
  real_type nrm2 = 0.0;
  for (sl_int i = 0; i < n; i++)
    nrm2 += (x.data())[i] * (x.data())[i];
  return sqrt(nrm2);
}

/* Only implements the simplest version of daxpy which is used in this code 
   here, with incrx and incry = 1 hardcoded. *y points at the results
   (overwrite) */
template <class real_type>
void daxpy(const sl_int n, const real alpha, const real_type x[],
	   const int incx, real_type y[], const int incy)
{
  real_type a = real_type(alpha);
  for (sl_int i = 0; i < n; i++)
    y[i] += a * x[i];
}

template <class real_type>
void daxpy(const sl_int n, const real alpha,
	   const boost::multi_array_ref<real_type, 3> &x, const int incx,
	   boost::multi_array_ref<real_type, 3> &y, const int incy)
{
  real_type a = real_type(alpha);
  for (sl_int i = 0; i < n; i++)
    (y.data())[i] += a * (x.data())[i];
}

template <class real_type>
void dcopy(const sl_int n, const real_type x[], const int incx, real_type y[],
	   const int incy)
{
  for (sl_int i = 0; i < n; i++)
    y[i] = x[i];
}

template <class real_type>
void dcopy(const sl_int n, const boost::multi_array_ref<real_type, 3> &x,
	   const int incx, boost::multi_array_ref<real_type, 3> &y, const int incy)
{
  for (sl_int i = 0; i < n; i++)
    (y.data())[i] = (x.data())[i];
}

#if defined(MKL_ILP64)

#include "mkl_cblas.h"

template <> inline
float ddot(const sl_int n, const float x[], const int incx,
	   const float y[], const int incy)
{
  return cblas_sdot(n, x, incx, y, incy);
}

template <> inline
double ddot(const sl_int n, const double x[], const int incx,
	    const double y[], const int incy)
{
  return cblas_ddot(n, x, incx, y, incy);
}

template <> inline
float ddot(const sl_int n, const boost::multi_array_ref<float, 3> &x,
	   const int incx, const boost::multi_array_ref<float, 3> &y,
	   const int incy)
{
  return cblas_sdot(n, x.data(), incx, y.data(), incy);
}

template <> inline
double ddot(const sl_int n, const boost::multi_array_ref<double, 3> &x,
	    const int incx, const boost::multi_array_ref<double, 3> &y,
	    const int incy)
{
  return cblas_ddot(n, x.data(), incx, y.data(), incy);
}

template <> inline
float dnrm2(const sl_int n, const float x[], const int incx)
{
  return cblas_snrm2(n, x, incx);
}

template <> inline
double dnrm2(const sl_int n, const double x[], const int incx)
{
  return cblas_dnrm2(n, x, incx);
}

template <> inline
float dnrm2(const sl_int n, const boost::multi_array_ref<float, 3> &x,
	    const int incx)
{
  return cblas_snrm2(n, x.data(), incx);
}

template <> inline
double dnrm2(const sl_int n, const boost::multi_array_ref<double, 3> &x,
	     const int incx)
{
  return cblas_dnrm2(n, x.data(), incx);
}

template <> inline
void daxpy(const sl_int n, const real alpha, const float x[],
	   const int incx, float y[], const int incy)
{
  cblas_saxpy(n, float(alpha), x, incx, y, incy);
}

template <> inline
void daxpy(const sl_int n, const real alpha, const double x[],
	   const int incx, double y[], const int incy)
{
  cblas_daxpy(n, double(alpha), x, incx, y, incy);
}

template <> inline
void daxpy(const sl_int n, const real alpha,
	   const boost::multi_array_ref<float, 3> &x, const int incx,
	   boost::multi_array_ref<float, 3> &y, const int incy)
{
  cblas_saxpy(n, float(alpha), x.data(), incx, y.data(), incy);
}

template <> inline
void daxpy(const sl_int n, const real alpha,
	   const boost::multi_array_ref<double, 3> &x, const int incx,
	   boost::multi_array_ref<double, 3> &y, const int incy)
{
  cblas_daxpy(n, double(alpha), x.data(), incx, y.data(), incy);
}

template <> inline
void dcopy(const sl_int n, const float x[], const int incx, float y[],
	   const int incy)
{
  cblas_scopy(n, x, incx, y, incy);
}

template <> inline
void dcopy(const sl_int n, const double x[], const int incx, double y[],
	   const int incy)
{
  cblas_dcopy(n, x, incx, y, incy);
}

template <> inline
void dcopy(const sl_int n, const boost::multi_array_ref<float, 3> &x,
	   const int incx, boost::multi_array_ref<float, 3> &y, const int incy)
{
  cblas_scopy(n, x.data(), incx, y.data(), incy);
}

template <> inline
void dcopy(const sl_int n, const boost::multi_array_ref<double, 3> &x,
	   const int incx, boost::multi_array_ref<double, 3> &y, const int incy)
{
  cblas_dcopy(n, x.data(), incx, y.data(), incy);
}

#endif // MKL

#endif // BLAS_WRAPPERS
