
#include <float.h>
#include "base_types.hpp"
#include "regularize.hpp"

static void D_func(const voxel_data &A, voxel_data &D1, voxel_data &D2,
		   voxel_data &D3, const int nx, const int ny, const int nz);
static void TV_main(const voxel_data &D1, const voxel_data &D2,
		    const voxel_data &D3, voxel_data &B, const int nx,
		    const int ny, const int nz);

#define EPS FLT_EPSILON //0.000001

inline static int sign(const double x)
{
  return (x > 0) - (x < 0);
}

void CCPi::tv_regularize(voxel_data &b, const voxel_data &a,
			 const int nx, const int ny, const int nz)
{
  voxel_data d1(boost::extents[nx][ny][nz]);
  voxel_data d2(boost::extents[nx][ny][nz]);
  voxel_data d3(boost::extents[nx][ny][nz]);
  D_func(a, d1, d2, d3, nx, ny, nz);
  TV_main(d1, d2, d3, b, nx, ny, nz);
}

// calculate differences 1,2,3
void D_func(const voxel_data &A, voxel_data &D1, voxel_data &D2,
	    voxel_data &D3, const int nx, const int ny, const int nz)
{
#pragma omp parallel for shared(A, D3) firstprivate(nx, ny, nz)
  for (int i = 0; i < nx; i++) {
    int i1 = i + 1;
    if (i1 >= nx)
      i1 = i - 1;
    int i2 = i - 1;
    if (i2 < 0)
      i2 = i + 1;
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
	  k1 = k-1;
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

	D1[i][j][k] = NOMx_1 / T1;
	D2[i][j][k] = NOMy_1 / T2;
	D3[i][j][k] = NOMz_1 / T3;
      }
    }
  }
}

// calculate divergence
void TV_main(const voxel_data &D1, const voxel_data &D2, const voxel_data &D3,
	     voxel_data &B, const int nx, const int ny, const int nz)
{
#pragma omp parallel for shared(B, D1, D2, D3) firstprivate(nx, ny, nz)
  for (int i = 0; i < nx; i++) {
    //int i1 = i + 1;
    //if (i1 >= nx)
    // i1 = i - 1;
    int i2 = i - 1;
    if (i2 < 0)
      i2 = i + 1;
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
	voxel_type dv1 = D1[i][j][k] - D1[i2][j][k];
	voxel_type dv2 = D2[i][j][k] - D2[i][j2][k];
	voxel_type dv3 = D3[i][j][k] - D3[i][j][k2];

	B[i][j][k] = dv1 + dv2 + dv3;
      }
    }
  }
}
