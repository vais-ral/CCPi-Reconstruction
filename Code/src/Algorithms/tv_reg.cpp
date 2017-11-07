
#include "base_types.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "blas.hpp"
#include "ui_calls.hpp"
#include "cgls.hpp"
#include "tv_reg.hpp"

// driver routine designed to initialise from CGLS
bool CCPi::tv_regularization::reconstruct(instrument *device,
					  voxel_data &voxels,
					  const real origin[3],
					  const real voxel_size[3])
{
  // Need copy of pixels before initial CGLS guess overwrites them
  pixel_data &p = device->get_pixel_data();
  int n_angles = device->get_num_angles();
  int n_h = device->get_num_h_pixels();
  int n_v = device->get_num_v_pixels();
  pixel_data b = pixel_data(boost::extents[n_angles][n_h][n_v]);
  copy(p, b, n_angles, n_h, n_v);
  cgls_3d cgls(5);
  bool ok = cgls.reconstruct(device, voxels, origin, voxel_size);
  if (ok)
    ok = reconstruct(device, b, voxels, origin, voxel_size);
  return ok;
}

// should be passed an initial guess in voxels, either from CGLS or
// some other method
bool CCPi::tv_regularization::reconstruct(instrument *device,
					  pixel_data &b, voxel_data &voxels,
					  const real origin[3],
					  const real voxel_size[3])
{
  const voxel_data::size_type *sz = voxels.shape();
  sl_int n_vox = sl_int(sz[0]) * sl_int(sz[1]) * sl_int(sz[2]);
  //voxel_type *const vox = voxels.data();

  // defaults - should be inputs
  real epsb_rel = 1e-4;
  int k_max = 10000;
  bool verbose = false;

  // Intitial settings of bmu and bL
  real bL = init_L;
  real bmu = init_mu;

  int numGrad = 0;
  int numBack = 0;
  int numFunc = 0;
  int numRest = 0;
  real_1d Lklist(k_max + 1);
  real_1d muklist(k_max + 1);

  // Todo - do we want to preserve the original voxels?
  //voxel_type *x = new voxel_type[n_vox];
  //for (sl_int i = 0; i < n_vox; i++)
  //x[i] = vox[i];
  real fxkp1 = 0.0;
  real hxkp1 = 0.0;
  real gxkp1 = 0.0;
  real_1d fxkp1l(k_max + 1);
  int k = 0;
			
  // Initialize vectors to hold tv and fidelity of the iterates
  // Todo - use this?
  bool ghxl = false;
  int tsize = 1;
  if (ghxl)
    tsize = k_max + 1;
  real_1d hxkp1l(tsize);
  real_1d gxkp1l(tsize);

  // Array to hold iterates in columns - Todo, use this?
  bool xl = false;
  real_1d xlist(1);
  if (xl and k_max * n_vox < 1e7)
    xlist.resize(n_vox * (k_max + 1));

  // setp up constraints - these are correct if constraint == 3
  real_1d c(1);
  real_1d d(1);
  if (constraint == 2) {
    // Todo - 3d voxel_data objects?
    c.resize(n_vox);
    d.resize(n_vox);
    for (sl_int i = 0; i < n_vox; i++)
      c[i] = 0.0;
    for (sl_int i = 0; i < n_vox; i++)
      d[i] = 1.0;
  } else {
    c[0] = 0.0;
    d[0] = 1.0;
  }

  std::list<int> rp;

  tvreg_core(voxels, fxkp1, hxkp1, gxkp1, fxkp1l, k, voxel_size, b, alpha,
	     tau, bL, bmu, epsb_rel, k_max, 3, sz[0], sz[1], sz[2], n_vox,
	     constraint, d, c, ghxl, xl, hxkp1l, gxkp1l, xlist, verbose,
	     numGrad, numBack, numFunc, numRest, Lklist, muklist, rp,
	     origin, device);

  // Todo - copy rp into rkList

  // If the iteration counter reaches the k_max, the algorithm did not converge
  if (k == k_max)
    report_error("Did not find a epsb_rel solution in k_max iterations.");
  return true;
}

bool CCPi::tv_regularization::supports_blocks() const
{
  return false;
}
