
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
  size_t id = get_global_id(0); 
  size_t jobid = get_global_id(1);
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
  size_t id = get_global_id(0); 
  size_t jobid = get_global_id(1);
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
  size_t id = get_global_id(0);
  size_t jobid = get_global_id(1);
  int xpos = ix * get_global_size(1) * xy_size + jobid * xy_size;
  // voxel ptr should be i/j column or needs an offset
  float vox = 0.0;
  int n = lengths[ix * get_global_size(1) + jobid];
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
  size_t id = get_global_id(0); 
  size_t jobid = get_global_id(1);
  int idx = id + id;
  int xpos = ix * get_global_size(1) * xy_size + jobid * xy_size;
  // voxel ptr should be i/j column or needs an offset
  float vox = 0.0;
  int n = lengths[ix * get_global_size(1) + jobid];
  for (int m = 0; m < n; m++) {
    // pixels needs m to be an a/h offset
    vox += (pixels[index[xpos + m] + id] + pixels[index[xpos + m] + id + 1])
      * l_xy[xpos + m];
  }
  voxels[offset + jobid * nz + id] += vox;
}

/*
__kernel void cone_xy_z_l(__global float *pixels, 
			  const __global float *voxels,
			  const __global float *alpha_xy,
			  const __global int *index,
			  const __global int *h,
			  const __global int *lengths,
			  const int nv, const int nz, const int start,
			  const int ah_size,
			  const float pzbz, const float inv_dz,
			  const int midp, const __global *delta_z,
			  const __global *inv_delz, const global *vox_z)
{
  // lower cone - 0 - < midp (size nv / 2)
  size_t id = get_global_id(0); // v
  size_t jobid = get_global_id(1);
  // pixel ptr should be ah column or needs an offset
  float pix = 0.0;
  int n = lengths[start + jobid];
  int pos = (start + jobid) * ah_size;

  // Todo - pass in from outside? 2D with pos?
  for (int l = 0; l < n; l++)
    alpha_inv[l] = alpha_xy[l] * inv_dz;

  float alpha_m0 = alpha_xy[pos];
  for (int m = 1; m < n; m++) {
    int k = int(pzbz + alpha_inv[pos + m - 1] * dz_ptr[id]);
    float alpha_m1 = alpha_xy[pos + m];
    float alpha_z = vz_ptr[k] * iz_ptr[id];
    float min_z = fmin(alpha_z, alpha_m1);
    pix += (voxels[index[pos + m] + k] * (min_z - alpha_m0)
	    + voxels[index[pos + m] + k - 1] * (alpha_m1 - min_z));
    alpha_m0 = alpha_m1;
  }
  pixels[h[start + jobid] * nv + id] += pix;  
}

__kernel void cone_xy_z_u(__global float *pixels, 
			  const __global float *voxels,
			  const __global float *alpha_xy,
			  const __global int *index,
			  const __global int *h,
			  const __global int *lengths,
			  const int nv, const int nz, const int start,
			  const int ah_size,
			  const float pzbz, const float inv_dz,
			  const int midp, const __global *delta_z,
			  const __global *inv_delz, const global *vox_z)
{
  // upper cone - midp - < nv (size nv / 2)
  size_t id = midp + get_global_id(0); // v
  size_t jobid = get_global_id(1);
  // pixel ptr should be ah column or needs an offset
  float pix = 0.0;
  int n = lengths[start + jobid];
  int pos = (start + jobid) * ah_size;

  // Todo - pass in from outside? 2D with pos?
  for (int l = 0; l < n; l++)
    alpha_inv[l] = alpha_xy[l] * inv_dz;

  float alpha_m0 = alpha_xy[pos];
  for (int m = 1; m < n; m++) {
    int k = int(pzbz + alpha_inv[pos + m - 1] * dz_ptr[id]);
    float alpha_m1 = alpha_xy[pos + m];
    float alpha_z = vz_ptr[k + 1] * iz_ptr[id];
    float min_z = fmin(alpha_z, alpha_m1);
    pix += (voxels[index[pos + m] + k] * (min_z - alpha_m0)
	    + voxels[index[pos + m] + k + 1] * (alpha_m1 - min_z));
    alpha_m0 = alpha_m1;
  }
  pixels[h[start + jobid] * nv + id] += pix;
}
*/

__kernel void cone_xy_z(__global float *pixels, 
			const __global float *voxels,
			const __global float *alpha_xy,
			const __global int *index,
			const __global int *h,
			const __global int *lengths,
			const int nv, const int nz, const int start,
			const int ah_size, const float pzbz,
			const int midp, const __global float *delta_z,
			const __global float *inv_delz,
			const global float *vox_z)
{
  size_t id = get_global_id(0); // v
  size_t jobid = get_global_id(1);
  // pixel ptr should be ah column or needs an offset
  float pix = 0.0;
  int n = lengths[start + jobid];
  int pos = (start + jobid) * ah_size;
  // hopefully this is -(1) + 0 or -0 + 1
  int vshift = -(id < midp) + (id >= midp);
  int zshift = (id >= midp);
  float del_z = delta_z[id];
  float inv_z = inv_delz[id];
  float alpha_m0 = alpha_xy[pos];
  for (int m = 1; m < n; m++) {
    int k = (int)(pzbz + alpha_m0 * del_z);
    float alpha_m1 = alpha_xy[pos + m];
    float alpha_z = vox_z[k + zshift] * inv_z;
    float min_z = fmin(alpha_z, alpha_m1);
    pix += (voxels[index[pos + m] + k] * (min_z - alpha_m0)
	    + voxels[index[pos + m] + k + vshift] * (alpha_m1 - min_z));
    alpha_m0 = alpha_m1;
  }
  pixels[h[start + jobid] * nv + id] += pix;
}

__kernel void cone_ah_z(const __global float *pixels, 
			__global float *voxels,
			const __global float *alpha_xy_0,
			const __global float *alpha_xy_1,
			const __global int *index,
			const __global int *h,
			const __global int *lengths,
			const int nv, const int nz, const int start,
			const int xy_size, const float pzbz,
			const int midp, const __global float *delta_z,
			const __global float *inv_delz,
			const global float *vox_z)
{
  size_t id = get_global_id(0); // v
  size_t jobid = get_global_id(1);
  // voxel ptr should be xy column or needs an offset
  int n = lengths[start + jobid];
  int pos = (start + jobid) * xy_size;
  // hopefully this is -(1) + 0 or -0 + 1
  int vshift = -(id < midp) + (id >= midp);
  int zshift = (id >= midp);
  float del_z = delta_z[id];
  float inv_z = inv_delz[id];
  for (int m = 0; m < n; m++) {
    const float alpha_m0 = alpha_xy_0[pos + m];
    int k = (int)(pzbz + alpha_m0 * del_z);
    float pix = pixels[index[pos + m] + id];
    const float alpha_m1 = alpha_xy_1[pos + m];
    float alpha_z = vox_z[k + zshift] * inv_z;
    float min_z = fmin(alpha_z, alpha_m1);
    // Todo - these need to be atomic - OpenCL 1.1 has no atomic_add(float)!
    voxels[h[start + jobid] * nz + k] += pix * (min_z - alpha_m0);
    voxels[h[start + jobid] * nz + k + vshift] += pix * (alpha_m1 - min_z);
  }
}
