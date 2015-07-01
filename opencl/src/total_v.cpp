
#include <float.h>
#include "base_types.hpp"
#include "regularize.hpp"

static void D1_func(const voxel_data &A, voxel_2d &D1, const int i,
		    const int nx, const int ny, const int nz);
static void D_func(const voxel_data &A, voxel_2d &D1, voxel_2d &D2,
		   voxel_2d &D3, const int i, const int nx, const int ny,
		   const int nz);
static void TV_main(const voxel_data &A, voxel_data &B, const int nx,
		    const int ny, const int nz);

#define EPS FLT_EPSILON //0.000001

inline static int sign(const double x)
{
  return (x > 0) - (x < 0);
}

void CCPi::tv_regularize(voxel_data &b, const voxel_data &a,
			 const int nx, const int ny, const int nz)
{
  TV_main(a, b, nx, ny, nz);
}

// calculate differences 1
void D1_func(const voxel_data &A, voxel_2d &D1, const int i,
	     const int nx, const int ny, const int nz)
{
  int i1 = i + 1;
  if (i1 >= nx)
    i1 = i - 1;
  int i2 = i - 1;
  if (i2 < 0)
    i2 = i + 1;
#pragma omp parallel for shared(A, D1) firstprivate(nx, ny, nz, i, i1, i2) schedule(dynamic)
  for (int j = 0; j < ny; j++) {
    int j1 = j + 1;
    if (j1 >= ny)
      j1 = j - 1;
    int j2 = j - 1;
    if (j2 < 0)
      j2 = j + 1;
    for (int k = 0; k < nz; k++) {
      // symmetric boundary conditions (Neuman)
      int k1 = k + 1;
      if (k1 >= nz)
	k1 = k - 1;
      int k2 = k - 1;
      if (k2 < 0)
	k2 = k + 1;

      // Forward-backward differences
      voxel_type NOMx_1 = A[i1][j][k] - A[i][j][k]; // x+
      voxel_type NOMy_1 = A[i][j1][k] - A[i][j][k]; // y+
      voxel_type NOMy_0 = A[i][j][k] - A[i][j2][k]; // y-
      voxel_type NOMz_1 = A[i][j][k1] - A[i][j][k]; // z+
      voxel_type NOMz_0 = A[i][j][k] - A[i][j][k2]; // z-

      voxel_type denom1x = NOMx_1 * NOMx_1;
      voxel_type denom2x = 0.5 * (sign(NOMy_1) + sign(NOMy_0))
	* std::min(std::abs(NOMy_1), std::abs(NOMy_0));
      denom2x = denom2x * denom2x;
      voxel_type denom3x = 0.5 * (sign(NOMz_1) + sign(NOMz_0))
	* std::min(std::abs(NOMz_1), std::abs(NOMz_0));
      denom3x = denom3x * denom3x;

      voxel_type T1 = std::sqrt(denom1x + denom2x + denom3x + EPS);

      D1[j][k] = NOMx_1 / T1;
    }
  }
}

// calculate differences 1,2,3
void D_func(const voxel_data &A, voxel_2d &D1, voxel_2d &D2,
	    voxel_2d &D3, const int i, const int nx, const int ny, const int nz)
{
  int i1 = i + 1;
  if (i1 >= nx)
    i1 = i - 1;
  int i2 = i - 1;
  if (i2 < 0)
    i2 = i + 1;
#pragma omp parallel for shared(A, D1, D2, D3) firstprivate(nx, ny, nz, i, i1, i2) schedule(dynamic)
  for (int j = 0; j < ny; j++) {
    int j1 = j + 1;
    if (j1 >= ny)
      j1 = j - 1;
    int j2 = j - 1;
    if (j2 < 0)
      j2 = j + 1;
    for (int k = 0; k < nz; k++) {
      // symmetric boundary conditions (Neuman)
      int k1 = k + 1;
      if (k1 >= nz)
	k1 = k - 1;
      int k2 = k - 1;
      if (k2 < 0)
	k2 = k + 1;

      // Forward-backward differences
      voxel_type NOMx_1 = A[i1][j][k] - A[i][j][k]; // x+
      voxel_type NOMy_1 = A[i][j1][k] - A[i][j][k]; // y+
      voxel_type NOMy_0 = A[i][j][k] - A[i][j2][k]; // y-
      voxel_type NOMx_0 = A[i][j][k] - A[i2][j][k]; // x-
      voxel_type NOMz_1 = A[i][j][k1] - A[i][j][k]; // z+
      voxel_type NOMz_0 = A[i][j][k] - A[i][j][k2]; // z-

      voxel_type denom1x = NOMx_1 * NOMx_1;
      voxel_type denom1y = NOMy_1 * NOMy_1;
      voxel_type denom1z = NOMz_1 * NOMz_1;

      voxel_type denom2x = 0.5 * (sign(NOMy_1) + sign(NOMy_0))
	* std::min(std::abs(NOMy_1), std::abs(NOMy_0));
      voxel_type denom2y = 0.5 * (sign(NOMx_1) + sign(NOMx_0))
	* std::min(std::abs(NOMx_1), std::abs(NOMx_0));
      //voxel_type denom2z = denom2y;

      denom2x = denom2x * denom2x;
      denom2y = denom2y * denom2y;
      //denom2z = denom2z * denom2z;
      voxel_type denom2z = denom2y;

      voxel_type denom3x = 0.5 * (sign(NOMz_1) + sign(NOMz_0))
	* std::min(std::abs(NOMz_1), std::abs(NOMz_0));
      //voxel_type denom3y = denom3x;
      //voxel_type denom3z = denom2x;

      denom3x = denom3x * denom3x;
      //denom3y = denom3y * denom3y;
      voxel_type denom3y = denom3x;
      //denom3z = denom3z * denom3z;
      voxel_type denom3z = denom2x;

      voxel_type T1 = std::sqrt(denom1x + denom2x + denom3x + EPS);
      voxel_type T2 = std::sqrt(denom1y + denom2y + denom3y + EPS);
      voxel_type T3 = std::sqrt(denom1z + denom2z + denom3z + EPS);

      D1[j][k] = NOMx_1 / T1;
      D2[j][k] = NOMy_1 / T2;
      D3[j][k] = NOMz_1 / T3;
    }
  }
}

// calculate divergence
void TV_main(const voxel_data &A, voxel_data &B,
	     const int nx, const int ny, const int nz)
{
  voxel_2d d1(boost::extents[ny][nz]);
  voxel_2d d2(boost::extents[ny][nz]);
  voxel_2d d3(boost::extents[ny][nz]);
  voxel_2d d1_i2(boost::extents[ny][nz]);
  // for start with i == 0, i2 == 1
  D1_func(A, d1_i2, 1, nx, ny, nz);
  for (int i = 0; i < nx; i++) {
    //int i1 = i + 1;
    //if (i1 >= nx)
    // i1 = i - 1;
    int i2 = i - 1;
    if (i2 < 0)
      i2 = i + 1;
    D_func(A, d1, d2, d3, i, nx, ny, nz);
#pragma omp parallel for shared(B, d1, d1_i2, d2, d3) firstprivate(nx, ny, nz, i, i2) schedule(dynamic)
    for (int j = 0; j < ny; j++) {
      //int j1 = j + 1;
      //if (j1 >= ny)
      //  j1 = j - 1;
      int j2 = j - 1;
      if (j2 < 0)
	j2 = j + 1;
      for (int k = 0; k < nz; k++) {
	// symmetric boundary conditions (Neuman)
	//int k1 = k + 1;
	//if (k1 >= nz)
	//  k1 = k-1;
	int k2 = k - 1;
	if (k2 < 0)
	  k2 = k + 1;

	// divergence components
	voxel_type dv1 = d1[j][k] - d1_i2[j][k];
	voxel_type dv2 = d2[j][k] - d2[j2][k];
	voxel_type dv3 = d3[j][k] - d3[j][k2];

	B[i][j][k] = dv1 + dv2 + dv3;
	// Could try be to smarter with double buffered pointers rather than
	// using a copy, for next i loop 'i2' is the 'i' value.
	d1_i2[j][k] = d1[j][k];
      }
    }
  }
}
