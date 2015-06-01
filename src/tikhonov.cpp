
#include "base_types.hpp"
#include "regularize.hpp"

void CCPi::tikhonov_regularize(voxel_data &b, const voxel_data &a,
			       const int nx, const int ny, const int nz)
{
#pragma omp parallel for shared(a, b) firstprivate(nx, ny, nz) schedule(dynamic)
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
      for (int k = 0; k < nz; k++)
	b[i][j][k] = 0.25 * (a[i1][j][k] + a[i2][j][k] + a[i][j1][k]
			     + a[i][j2][k]) - a[i][j][k];
    }
  }
}
