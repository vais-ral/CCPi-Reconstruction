
#include "base_types.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "ui_calls.hpp"

// driver routine designed to initialise from CGLS
bool CCPi::tv_regularization::reconstruct(const instrument *device,
					  voxel_data &voxels,
					  const real origin[3],
					  const real voxel_size[3])
{
  // Need copy of pixels before initial CGLS guess overwrites them
  sl_int n_pixels = device->get_data_size();
  pixel_type *const p = device->get_pixel_data();
  pixel_type *b = new pixel_type[n_pixels];
  for (sl_int i = 0; i < n_pixels; i++)
    b[i] = p[i];
  cgls_base cgls(5);
  bool ok = cgls.reconstruct(device, voxels, origin, voxel_size);
  if (ok)
    ok = reconstruct(device, b, voxels, origin, voxel_size);
  delete [] p;
  return ok;
}

// should be passed an initial guess in voxels, either from CGLS or
// some other method
bool CCPi::tv_regularization::reconstruct(const instrument *device,
					  pixel_type *b, voxel_data &voxels,
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
  std::vector<real> Lklist(k_max + 1);
  std::vector<real> muklist(k_max + 1);

  // Todo - do we want to preserve the original voxels?
  //voxel_type *x = new voxel_type[n_vox];
  //for (sl_int i = 0; i < n_vox; i++)
  //x[i] = vox[i];
  real fxkp1 = 0.0;
  real hxkp1 = 0.0;
  real gxkp1 = 0.0;
  std::vector<real> fxkp1l(k_max + 1);
  int k = 0;
			
  // Initialize vectors to hold tv and fidelity of the iterates
  // Todo - use this?
  bool ghxl = false;
  int tsize = 1;
  if (ghxl)
    tsize = k_max + 1;
  std::vector<real> hxkp1l(tsize);
  std::vector<real> gxkp1l(tsize);

  // Array to hold iterates in columns - Todo, use this?
  bool xl = false;
  std::vector<real> xlist(1);
  if (xl and k_max * n_vox < 1e7)
    xlist.resize(n_vox * (k_max + 1));

  // setp up constraints - these are correct if constraint == 3
  std::vector<real> c(1);
  std::vector<real> d(1);
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
