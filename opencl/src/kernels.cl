
// OpenCL kernels - only use 2 dims for global workspace - should ids be local?
// should l_xy/index be __constant?

__kernel void parallel_xy_z1(__global float *pixels, const int offset,
			     const __global float *voxels,
			     const __global float *l_xy,
			     const __global int *index, const int n,
			     const int nv, const int nz)
{
  // ? Index in work set to vectorize z - nz/nv total
  size_t id = get_global_id(1) * get_global_size(0) + get_global_id(0); 
  // pixel ptr should be ah column or needs an offset
  float pix = 0.0;
  for (int m = 0; m < n; m++) {
    // voxels needs m to be an i/j offset
    pix += voxels[index[m] + id] * l_xy[m];
  }
  pixels[offset + id] += pix;
}

__kernel void parallel_xy_z2(__global float *pixels, const int offset,
			     const __global float *voxels,
			     const __global float *l_xy,
			     const __global int *index, const int n,
			     const int nv, const int nz)
{
  // ? Index in work set to vectorize z - nz total!
  size_t id = get_global_id(1) * get_global_size(0) + get_global_id(0); 
  // pixel ptr should be ah column or needs an offset
  int idx = id + id;
  float pix1 = 0.0;
  float pix2 = 0.0;
  for (int m = 0; m < n; m++) {
    // voxels needs m to be an i/j offset
    pix1 += voxels[index[m] + id] * l_xy[m];
    pix2 += voxels[index[m] + id] * l_xy[m];
  }
  pixels[offset + idx] += pix1;
  pixels[offset + idx + 1] += pix2;
}

__kernel void parallel_ah_z1(const __global float *pixels,
			     __global float *voxels, const int offset,
			     const __global float *l_xy,
			     const __global int *index, const int n,
			     const int nv, const int nz)
{
  // ? Index in work set to vectorize z - nz total!
  size_t id = get_global_id(1) * get_global_size(0) + get_global_id(0); 
  // voxel ptr should be i/j column or needs an offset
  float vox = 0.0;
  for (int m = 0; m < n; m++) {
    // pixels needs m to be an a/h offset
    vox += pixels[index[m] + id] * l_xy[m];
  }
  voxels[offset + id] += vox;
}

__kernel void parallel_ah_z2(const __global float *pixels,
			     __global float *voxels, const int offset,
			     const __global float *l_xy,
			     const __global int *index, const int n,
			     const int nv, const int nz)
{
  // ? Index in work set to vectorize z - nz total!
  size_t id = get_global_id(1) * get_global_size(0) + get_global_id(0); 
  int idx = id + id;
  // voxel ptr should be i/j column or needs an offset
  float vox = 0.0;
  for (int m = 0; m < n; m++) {
    // pixels needs m to be an a/h offset
    vox += (pixels[index[m] + idx] + pixels[index[m] + idx + 1]) * l_xy[m];
  }
  voxels[offset + id] += vox;
}
