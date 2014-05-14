
#ifndef BLAS_WRAPPERS
#define BLAS_WRAPPERS

template <class real_type>
inline void init_data(boost::multi_array_ref<real_type, 3> &x,
		      const sl_int nx, const sl_int ny, const sl_int nz)
{
  sl_int n = nx * ny * nz;
#pragma omp parallel for shared(x) firstprivate(n) schedule(dynamic)
  for (sl_int i = 0; i < n; i++)
    (x.data())[i] = 0.0;
}

template <class real_type>
inline void scal_y(const real beta, boost::multi_array_ref<real_type, 3> &y,
		   const sl_int nx, const sl_int ny, const sl_int nz)
{
  // y = beta * y
  sl_int n = nx * ny * nz;
  real_type b = real_type(beta);
#pragma omp parallel for shared(y) firstprivate(n, b) schedule(dynamic)
  for (sl_int i = 0; i < n; i++)
    (y.data())[i] = b * (y.data())[i];
}

inline voxel_type norm_voxels(const voxel_data &v, const sl_int nx,
			      const sl_int ny, const sl_int nz)
{
  voxel_type norm = 0.0;
#pragma omp parallel for reduction(+:norm) shared(v) firstprivate(nx, ny, nz) schedule(dynamic)
  for (sl_int i = 0; i < nz; i++) {
    voxel_type n1 = 0.0;
    for (sl_int j = 0; j < ny; j++) {
      voxel_type n2 = 0.0;
      for (sl_int k = 0; k < nx; k++)
	n2 += v[k][j][i] * v[k][j][i];
      n1 += n2;
    }
    norm += n1;
  }
  return norm;
}

inline void norm_voxels(const voxel_data &v, const sl_int nx, const sl_int ny,
			const sl_int nz, voxel_1d &norm)
{
#pragma omp parallel for shared(v) firstprivate(nx, ny, nz) schedule(dynamic)
  for (sl_int i = 0; i < nz; i++) {
    voxel_type n1 = 0.0;
    for (sl_int j = 0; j < ny; j++) {
      voxel_type n2 = 0.0;
      for (sl_int k = 0; k < nx; k++)
	n2 += v[k][j][i] * v[k][j][i];
      n1 += n2;
    }
    norm[i] = n1;
  }
}

inline pixel_type norm_pixels(const pixel_data &p, const sl_int na,
			      const sl_int nv, const sl_int nh)
{
  pixel_type norm = 0.0;
#pragma omp parallel for reduction(+:norm) shared(p) firstprivate(na, nv, nh) schedule(dynamic)
  for (sl_int i = 0; i < na; i++) {
    pixel_type n1 = 0.0;
    for (sl_int j = 0; j < nv; j++) {
      pixel_type n2 = 0.0;
      for (sl_int k = 0; k < nh; k++)
	n2 += p[i][j][k] * p[i][j][k];
      n1 += n2;
    }
    norm += n1;
  }
  return norm;
}

inline void norm_pixels(const pixel_data &p, const sl_int na, const sl_int nv,
			const sl_int nh, pixel_1d &norm)
{
#pragma omp parallel for shared(p) firstprivate(na, nv, nh) schedule(dynamic)
  for (sl_int j = 0; j < nv; j++) {
    pixel_type n1 = 0.0;
    for (sl_int i = 0; i < na; i++) {
      pixel_type n2 = 0.0;
      for (sl_int k = 0; k < nh; k++)
	n2 += p[i][j][k] * p[i][j][k];
      n1 += n2;
    }
    norm[j] = n1;
  }
}

template <class real_type>
inline void sum_axpy(const real alpha,
		     const boost::multi_array_ref<real_type, 3> &x,
		     boost::multi_array_ref<real_type, 3> &y, const sl_int nx,
		     const sl_int ny, const sl_int nz)
{
  // y += alpha * x
  sl_int n = nx * ny * nz;
  real_type a = real_type(alpha);
#pragma omp parallel for shared(x, y) firstprivate(n, a) schedule(dynamic)
  for (sl_int i = 0; i < n; i++)
    (y.data())[i] += a * (x.data())[i];
}

template <class real_type>
inline void sum_axpy(const pixel_1d &alpha,
		     const boost::multi_array_ref<real_type, 3> &x,
		     boost::multi_array_ref<real_type, 3> &y, const sl_int nx,
		     const sl_int ny, const sl_int nz)
{
  // y += alpha * x
#pragma omp parallel for shared(x, y, alpha) firstprivate(nx, ny, nz) schedule(dynamic)
  for (sl_int i = 0; i < nz; i++) {
    voxel_type a = alpha[i];
    for (sl_int j = 0; j < ny; j++)
      for (sl_int k = 0; k < nx; k++)
	y[k][j][i] += a * x[k][j][i];
  }
}

template <class real_type>
inline void sub_axpy(const pixel_1d &alpha,
		     const boost::multi_array_ref<real_type, 3> &x,
		     boost::multi_array_ref<real_type, 3> &y, const sl_int nx,
		     const sl_int ny, const sl_int nz, const int ppv)
{
  // y += alpha * x
#pragma omp parallel for shared(x, y, alpha) firstprivate(nx, ny, nz, ppv) schedule(dynamic)
  for (sl_int i = 0; i < nz; i++) {
    voxel_type a = alpha[i / ppv];
    for (sl_int j = 0; j < ny; j++)
      for (sl_int k = 0; k < nx; k++)
	y[k][j][i] -= a * x[k][j][i];
  }
}

template <class real_type>
inline void scal_xby(const boost::multi_array_ref<real_type, 3> &x,
		     const real beta, boost::multi_array_ref<real_type, 3> &y,
		     const sl_int nx, const sl_int ny, const sl_int nz)
{
  // y = x + beta * y
  sl_int n = nx * ny * nz;
  real_type b = real_type(beta);
#pragma omp parallel for shared(x, y) firstprivate(n, b) schedule(dynamic)
  for (sl_int i = 0; i < n; i++)
    (y.data())[i] = b * (y.data())[i] + (x.data())[i];
}

template <class real_type>
inline void scal_xby(const boost::multi_array_ref<real_type, 3> &x,
		     const voxel_1d &beta,
		     boost::multi_array_ref<real_type, 3> &y,
		     const sl_int nx, const sl_int ny, const sl_int nz)
{
  // y = x + beta * y
#pragma omp parallel for shared(x, y, beta) firstprivate(nx, ny, nz) schedule(dynamic)
  for (sl_int i = 0; i < nz; i++) {
    voxel_type b = beta[i];
    for (sl_int j = 0; j < ny; j++)
      for (sl_int k = 0; k < nx; k++)
	y[k][j][i] = b * y[k][j][i] + x[k][j][i];
  }
}

template <class real_type>
inline void sum_xbyz(const boost::multi_array_ref<real_type, 3> &x,
		     const real beta,
		     const boost::multi_array_ref<real_type, 3> &y,
		     boost::multi_array_ref<real_type, 3> &z,
		     const sl_int nx, const sl_int ny, const sl_int nz)
{
  // z = x + b * y
  sl_int n = nx * ny * nz;
  real_type b = real_type(beta);
#pragma omp parallel for shared(x, y) firstprivate(n, b) schedule(dynamic)
  for (sl_int i = 0; i < n; i++)
    (z.data())[i] = (x.data())[i] + b * (y.data())[i];
}

template <class real_type>
inline void diff_xyz(const boost::multi_array_ref<real_type, 3> &x,
		     const boost::multi_array_ref<real_type, 3> &y,
		     boost::multi_array_ref<real_type, 3> &z,
		     const sl_int nx, const sl_int ny, const sl_int nz)
{
  sum_xbyz(x, -1.0, y, z, nx, ny, nz);
  /*
  // z = x - y
  sl_int n = nx * ny * nz;
#pragma omp parallel for shared(x, y) firstprivate(n) schedule(dynamic)
  for (sl_int i = 0; i < n; i++)
    (z.data())[i] = (x.data())[i] - (y.data())[i];
  */
}

inline voxel_type dot_prod(const voxel_data &x, const voxel_data &y,
			   const sl_int nx, const sl_int ny, const sl_int nz)
{
  voxel_type norm = 0.0;
#pragma omp parallel for reduction(+:norm) shared(x, y) firstprivate(nx, ny, nz) schedule(dynamic)
  for (sl_int i = 0; i < nz; i++) {
    voxel_type n1 = 0.0;
    for (sl_int j = 0; j < ny; j++) {
      voxel_type n2 = 0.0;
      for (sl_int k = 0; k < nx; k++)
	n2 += x[k][j][i] * y[k][j][i];
      n1 += n2;
    }
    norm += n1;
  }
  return norm;
}

template <class real_type>
inline void copy(const boost::multi_array_ref<real_type, 3> &x,
		 boost::multi_array_ref<real_type, 3> &y,
		 const sl_int nx, const sl_int ny, const sl_int nz)
{
  // z = x - y
  sl_int n = nx * ny * nz;
#pragma omp parallel for shared(x, y) firstprivate(n) schedule(dynamic)
  for (sl_int i = 0; i < n; i++)
    (y.data())[i] = (x.data())[i];
}

#if defined(MKL_ILP64)

#include "mkl_cblas.h"

template <>
inline void scal_y(const real beta, boost::multi_array_ref<float, 3> &y,
		   const sl_int nx, const sl_int ny, const sl_int nz)
{
  sl_int n = nx * ny * nz;
  cblas_sscal(n, float(beta), y.data(), 1);
}

template <>
inline void scal_y(const real beta, boost::multi_array_ref<double, 3> &y,
		   const sl_int nx, const sl_int ny, const sl_int nz)
{
  sl_int n = nx * ny * nz;
  cblas_dscal(n, double(beta), y.data(), 1);
}

template <>
inline void sum_axpy(const real alpha,
		     const boost::multi_array_ref<float, 3> &x,
		     boost::multi_array_ref<float, 3> &y, const sl_int nx,
		     const sl_int ny, const sl_int nz)
{
  // y += alpha * x
  sl_int n = nx * ny * nz;
  cblas_saxpy(n, float(alpha), x.data(), 1, y.data(), 1);
}

template <>
inline void sum_axpy(const real alpha,
		     const boost::multi_array_ref<double, 3> &x,
		     boost::multi_array_ref<double, 3> &y, const sl_int nx,
		     const sl_int ny, const sl_int nz)
{
  // y += alpha * x
  sl_int n = nx * ny * nz;
  cblas_daxpy(n, double(alpha), x.data(), 1, y.data(), 1);
}

template <> inline
void copy(const boost::multi_array_ref<float, 3> &x,
	  boost::multi_array_ref<float, 3> &y, const sl_int nx,
	  const sl_int ny, const sl_int nz)
{
  cblas_scopy(nx * ny * nz, x.data(), 1, y.data(), 1);
}

template <> inline
void copy(const boost::multi_array_ref<double, 3> &x,
	  boost::multi_array_ref<double, 3> &y, const sl_int nx,
	  const sl_int ny, const sl_int nz)
{
  cblas_dcopy(nx * ny * nz, x.data(), 1, y.data(), 1);
}

#endif // MKL

#endif // BLAS_WRAPPERS
