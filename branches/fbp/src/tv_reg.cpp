
#include "base_types.hpp"
#include "fbp.hpp"
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
  voxel_type *const vox = voxels.data();

  // defaults - should be inputs
  real epsb_rel = 1e-4;
  int k_max = 10000;
  bool verbose = false;

  // Intitial settings of bmu and bL
  real bL = init_L;
  real bmu = init_mu;

  real numGrad = 0.0;
  real numBack = 0.0;
  real numFunc = 0.0;
  real numRest = 0.0;
  real *Lklist = new real[k_max + 1];
  real *muklist = new real[k_max + 1];

  // Todo - do we want to preserve the original voxels?
  //voxel_type *x = new voxel_type[n_vox];
  //for (sl_int i = 0; i < n_vox; i++)
  //x[i] = vox[i];
  real fxkp1 = 0.0;
  real hxkp1 = 0.0;
  real gxkp1 = 0.0;
  real *fxkp1l = new real[k_max + 1];
  int k = 0;
			
  // Initialize vectors to hold tv and fidelity of the iterates
  // Todo - use this?
  bool ghxl = false;
  real *hxkp1l = 0;
  real *gxkp1l = 0;
  if (ghxl) {
    hxkp1l = new real[k_max + 1];
    gxkp1l = new real[k_max + 1];
  }

  // Array to hold iterates in columns - Todo, use this?
  bool xl = false;
  real *xlist = 0;
  if (xl and k_max * n_vox < 1e7)
    xlist = new real[n_vox * (k_max + 1)];

  // setp up constraints - these are correct if constraint == 3
  real cbase = 0.0;
  real dbase = 1.0;
  real *c = &cbase;
  real *d = &dbase;
  if (constraint == 2) {
    c = new real[n_vox];
    d = new real[n_vox];
    for (sl_int i = 0; i < n_vox; i++)
      c[i] = 0.0;
    for (sl_int i = 0; i < n_vox; i++)
      d[i] = 1.0;
  }

  std::list<int> rp;

  tvreg_core(vox, &fxkp1, &hxkp1, &gxkp1, fxkp1l, &k, voxel_size, b, alpha,
	     tau, bL, bmu, epsb_rel, k_max, 3, sz[0], sz[1], sz[2], n_vox,
	     constraint, d, c, ghxl, xl, hxkp1l, gxkp1l, xlist, verbose,
	     &numGrad, &numBack, &numFunc, &numRest, Lklist, muklist, rp,
	     origin, device);

  if (constraint == 2) {
    delete [] c;
    delete [] d;
  }

  // Todo - use xlist etc and other values
  // Truncate excess zeros in fxkp1l, xlist, gxkp1l and hxkp1l
  if (ghxl) {
    delete [] gxkp1l;
    delete [] hxkp1l;
  }
  if (xlist != 0)
    delete [] xlist;

  delete [] fxkp1l;
  //delete [] x;
  delete [] muklist;
  delete [] Lklist;

  // Todo - copy rp into rkList

  // If the iteration counter reaches the k_max, the algorithm did not converge
  if (k == k_max)
    report_error("Did not find a epsb_rel solution in k_max iterations.");
  return true;
}
