
// OpenCL kernels - only use 2 dims for global workspace - should ids be local?
// should l_xy/index be __constant?

__kernel void parallel_xy_z1(__global float *pixels,
			     const __global float *voxels,
			     const __global float *l_xy,
			     const __global int *index,
			     const __global int *h,
			     const __global int *lengths,
			     const int nv, const int nz, const int start,
			     const int ah_size)
{
  // ? Index in work set to vectorize z - nz/nv total
  size_t id = get_global_id(1) * get_global_size(0) + get_global_id(0); 
  size_t jobid = get_global_id(2);
  // pixel ptr should be ah column or needs an offset
  float pix = 0.0;
  int n = lengths[start + jobid];
  int pos = (start + jobid) * ah_size;
  for (int m = 0; m < n; m++) {
    // voxels needs m to be an i/j offset
    pix += voxels[index[pos + m] + id] * l_xy[pos + m];
  }
  pixels[h[start + jobid] * nv + id] += pix;
}

__kernel void parallel_xy_z2(__global float *pixels,
			     const __global float *voxels,
			     const __global float *l_xy,
			     const __global int *index,
			     const __global int *h,
			     const __global int *lengths,
			     const int nv, const int nz, const int start,
			     const int ah_size)
{
  // ? Index in work set to vectorize z - nz total!
  size_t id = get_global_id(1) * get_global_size(0) + get_global_id(0); 
  size_t jobid = get_global_id(2);
  // pixel ptr should be ah column or needs an offset
  int idx = id + id;
  int n = lengths[start + jobid];
  int pos = (start + jobid) * ah_size;
  float pix1 = 0.0;
  float pix2 = 0.0;
  for (int m = 0; m < n; m++) {
    // voxels needs m to be an i/j offset
    pix1 += voxels[index[pos + m] + id] * l_xy[pos + m];
    pix2 += voxels[index[pos + m] + id] * l_xy[pos + m];
  }
  pixels[h[start + jobid] * nv + idx] += pix1;
  pixels[h[start + jobid] * nv + idx + 1] += pix2;
}

__kernel void parallel_ah_z1(const __global float *pixels,
			     __global float *voxels, const int offset,
			     const __global float *l_xy,
			     const __global int *index,
			     const __global int *lengths,
			     const int nv, const int nz, const int xy_size,
			     const int ix)
{
  // ? Index in work set to vectorize z - nz total!
  size_t id = get_global_id(1) * get_global_size(0) + get_global_id(0);
  size_t jobid = get_global_id(2);
  int xpos = ix * get_global_size(2) * xy_size + jobid * xy_size;
  // voxel ptr should be i/j column or needs an offset
  float vox = 0.0;
  int n = lengths[ix * get_global_size(2) + jobid];
  for (int m = 0; m < n; m++) {
    // pixels needs m to be an a/h offset
    vox += pixels[index[xpos + m] + id] * l_xy[xpos + m];
  }
  voxels[offset + jobid * nz + id] += vox;
}

__kernel void parallel_ah_z2(const __global float *pixels,
			     __global float *voxels, const int offset,
			     const __global float *l_xy,
			     const __global int *index,
			     const __global int *lengths,
			     const int nv, const int nz, const int xy_size,
			     const int ix)
{
  // ? Index in work set to vectorize z - nz total!
  size_t id = get_global_id(1) * get_global_size(0) + get_global_id(0); 
  size_t jobid = get_global_id(2);
  int idx = id + id;
  int xpos = ix * get_global_size(2) * xy_size + jobid * xy_size;
  // voxel ptr should be i/j column or needs an offset
  float vox = 0.0;
  int n = lengths[ix * get_global_size(2) + jobid];
  for (int m = 0; m < n; m++) {
    // pixels needs m to be an a/h offset
    vox += (pixels[index[xpos + m] + id] + pixels[index[xpos + m] + id + 1])
      * l_xy[xpos + m];
  }
  voxels[offset + jobid * nz + id] += vox;
}
