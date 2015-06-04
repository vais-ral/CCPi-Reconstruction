
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
      const voxel_type *a_i1 = assume_aligned(&(a[i1][j][0]), voxel_type);
      const voxel_type *a_i2 = assume_aligned(&(a[i2][j][0]), voxel_type);
      const voxel_type *a_j1 = assume_aligned(&(a[i][j1][0]), voxel_type);
      const voxel_type *a_j2 = assume_aligned(&(a[i][j2][0]), voxel_type);
      const voxel_type *a_ij = assume_aligned(&(a[i][j][0]), voxel_type);
      voxel_type *b_ij = assume_aligned(&(b[i][j][0]), voxel_type);
      for (int k = 0; k < nz; k++)
	b_ij[k] = 0.25 * (a_i1[k] + a_i2[k] + a_j1[k] + a_j2[k]) - a_ij[k];
    }
  }
}
