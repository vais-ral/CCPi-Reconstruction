
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#include <float.h>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // types
#include "blas.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "ui_calls.hpp"
#include "Algorithms/tv_reg.hpp"

/* Settings which makes the user do a CTRL-C break out of the loop*/
#if defined(_WIN32) || defined(__WIN32__)

#  ifdef LIBUT
#    define STOPMARK utIsInterruptPending()
#    define INITBREAK ;
bool utIsInterruptPending(void);
#  else
#    define STOPMARK false
#    define INITBREAK ;
#  endif // LIBUT

#else

#include <signal.h>
#define INITBREAK   sigint_cont = 1;  (void) signal(SIGINT , ex_sigint);
#define STOPMARK sigint_cont==0
int sigint_cont = 1;
void ex_sigint(int sig) {
	sigint_cont = 0;
}
#endif

static void P(voxel_data &y, const int ctype, const real_1dr &d,
	      const real_1dr &c, const voxel_data::size_type sz[]);
static real DTD(voxel_data &x, voxel_data &Nablafx, std::vector<real> &uijl,
		const real tau, const int Ddim, const int Dm, const int Dn,
		const int Dl, const voxel_data::size_type sz[]);

void CCPi::tv_regularization::tvreg_core(voxel_data &xkp1, real &fxkp1,
					 real &hxkp1, real &gxkp1,
					 real_1dr &fxkp1l,
					 int &kend, const real voxel_size[],
					 const pixel_data &b, const real alpha,
					 real tau, real bL, real bmu,
					 real epsb_rel, int k_max,
					 const int Ddim, const int Dm,
					 const int Dn, const int Dl,
					 const sl_int prodDims, const int ctype,
					 real_1dr &d,
					 real_1dr &c, const bool ghxl,
					 const bool xl,
					 real_1dr &hxkp1l,
					 real_1dr &gxkp1l,
					 real_1dr &xlist,
					 const bool verbose, int &numGrad,
					 int &numBack, int &numFunc,
					 int &numRest,
					 real_1dr &Lklist,
					 real_1dr &muklist,
					 std::list<int> &rp,
					 const real grid_offset[],
					 instrument *device)
{
  // FLT_EPSILON? or 1e-14
  const real one_eps = 1.0 + 1e-14; //1.0 + FLT_EPSILON;
  real q;
  real thetakp1;
  real betak;
  real nGt;
  real t;
  real s_L = 1.3;
  real s_mu = 0.7;
  real Lm1;
  real nGtm1;
  real gamma0 = 0.0;
  int n_angles = device->get_num_angles();
  int n_v = device->get_num_v_pixels();
  int n_h = device->get_num_h_pixels();
  //sl_int n_rays = sl_int(n_angles) * sl_int(n_v) * sl_int(n_h);

  INITBREAK

  const voxel_data::size_type *sz = xkp1.shape();
  //sl_int n_vox = sl_int(sz[0]) * sl_int(sz[1]) * sl_int(sz[2]);
  voxel_3d Nablafyk(boost::extents[sz[0]][sz[1]][sz[2]],
		    boost::c_storage_order());
  voxel_3d Nablafxkp1(boost::extents[sz[0]][sz[1]][sz[2]],
		      boost::c_storage_order());
  voxel_3d yk(boost::extents[sz[0]][sz[1]][sz[2]],
	      boost::c_storage_order());
  voxel_3d xk(boost::extents[sz[0]][sz[1]][sz[2]],
	      boost::c_storage_order());
  voxel_3d temp(boost::extents[sz[0]][sz[1]][sz[2]],
		boost::c_storage_order());
  sl_int nx = sl_int(sz[0]);
  sl_int ny = sl_int(sz[1]);
  sl_int nz = sl_int(sz[2]);

  /*temp vectors */
  voxel_3d tv(boost::extents[sz[0]][sz[1]][sz[2]],
	      boost::c_storage_order());
  pixel_3d tv2(boost::extents[n_angles][n_h][n_v]);
  std::vector<real> uijl(Ddim);

  add_output("TVreg initialising...");
  send_output();

  /* INITIALIZE */
  numGrad = 0;
  numBack = 0;
  numFunc = 0;
  numRest = 0;

 restart:
  real cumprod = 1.0;

  /* Project solution onto feasible space */
  P(xkp1, ctype, d, c, sz);

  real thetak = std::sqrt(bmu / bL);

  /* Include a backtracking to initialize everything */
  /* y_k = x_k+1 */
  copy(xkp1, yk, nx, ny, nz);

  /* Calculate the gradient in yk=xk+1 */
  numGrad++;

  /* alpha*T_tau(y_k) (Nablafyk is returned as gradient) */
  real fyk = alpha * DTD(yk, Nablafyk, uijl, tau, Ddim, Dm, Dn, Dl, sz);

  /* For each voxel */
  scal_y(alpha, Nablafyk, nx, ny, nz);
  /*-----------------Forward projection--------------------------------*/
  init_data(tv2, n_angles, n_h, n_v);
  /*tv2 = A*yk */
  device->forward_project(tv2, yk, grid_offset, voxel_size, Dm, Dn, Dl);
  /* tv2 = tv2 - b */
  sum_axpy(real(-1.0), b, tv2, n_angles, n_h, n_v);

  numFunc++;
  /* fyk + 0.5*||A*y_k - b||^2 */
  fyk += real(0.5) * norm_pixels(tv2, n_angles, n_v, n_h);

  init_data(tv, nx, ny, nz);
  /*------------------Backward projection------------------------------*/
  init_data(temp, nx, ny, nz);
  device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
  /* For each voxel */
  sum_axpy(1.0, temp, Nablafyk, nx, ny, nz);

  /* Take the projected step from yk to xkp1 */
  /* bL is original setting of L_k */
  t = - 1 / bL;
  /* x_k+1 = y_k - t*Nablaf(y_k) */
  sum_xbyz(yk, t, Nablafyk, xkp1, nx, ny, nz);

  P(xkp1, ctype, d, c, sz);

  /* Backtracking on Lipschitz parameter. */
  diff_xyz(xkp1, yk, tv, nx, ny, nz);

  /* alpha*T_tau(xkp1+1) (Nablafyk is returned as gradient) */
  hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl, sz);

  /*-----------------Forward projection--------------------------------*/
  init_data(tv2, n_angles, n_h, n_v);
  /*tv2 = A*x_k+1 */
  device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
  /* tv2 = (A*x_k+1 - b) */
  sum_axpy(-1.0, b, tv2, n_angles, n_h, n_v);

  /* 0.5*||A*x_k+1 - b||^2 */
  gxkp1 = real(0.5) * norm_pixels(tv2, n_angles, n_v, n_h);

  numFunc++;
  // f(x_k+1) = h(x_k+1) + g(x_k+1) = alpha*T_tau(x_k+1) + 0.5*||A*x_k+1 - b||^2
  fxkp1 = hxkp1 + gxkp1;

  while (fxkp1 / one_eps > fyk + dot_prod(Nablafyk, tv, nx, ny, nz)
	 + (bL / 2) * norm_voxels(tv, nx, ny, nz)) {
    numBack++;
    bL *= s_L;

    /* Take the projected step from yk to xkp1 */
    /* t = - L^-1 */
    t = - 1 / bL;
    /* x_k+1 = y_k - t*Nablafyk */
    sum_xbyz(yk, t, Nablafyk, xkp1, nx, ny, nz);
    /* project to within bounds (x_k+1) */
    P(xkp1, ctype, d, c, sz);

    /* Backtracking on Lipschitz parameter. */
    diff_xyz(xkp1, yk, tv, nx, ny, nz);

    /* alpha*T_tau(x_k+1) */
    hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl, sz);

    /*-----------------Forward projection--------------------------------*/
    init_data(tv2, n_angles, n_h, n_v);
    /*tv2 = A*x_k+1 */
    device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
    /* tv2 = (A*x_k+1 - b) */
    sum_axpy(-1.0, b, tv2, n_angles, n_h, n_v);

    /* 0.5*||A*x_k+1 - b||^2 */
    gxkp1 = real(0.5) * norm_pixels(tv2, n_angles, n_v, n_h);

    numFunc++;
    /* f(x_k+1) = h(x_k+1) + g(x_k+1)
       = alpha*T_tau(x_k+1) + 0.5*||A*x_k+1 - b||^2 */
    fxkp1 = hxkp1 + gxkp1;
  }

  /* Calculate initial gradient map */
  if (ctype == 1) { // No constraints - just real numbers
    nGt = std::sqrt(norm_voxels(Nablafxkp1, nx, ny, nz));
  } else { /* Positivity constraint */
    t = - 1 / bL;
    sum_xbyz(xkp1, t, Nablafxkp1, tv, nx, ny, nz);
    P(tv, ctype, d, c, sz);

    scal_xby(xkp1, real(-1.0), tv, nx, ny, nz);

    nGt = bL * std::sqrt(norm_voxels(tv, nx, ny, nz));
  }

  /* save the initial parameter */
  Lm1 = bL;
  nGtm1 = nGt;
  copy(xkp1, xk, nx, ny, nz);
  copy(xkp1, yk, nx, ny, nz);

  /* LOOP */
  bool stop = false; /*Flag for when to break the for-loop*/

  /*Calculate fxk */
  real fxk = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,Ddim,Dm,Dn,Dl, sz);

  /*-----------------Forward projection--------------------------------*/
  init_data(tv2, n_angles, n_h, n_v);
  device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
  sum_axpy(real(-1.0), b, tv2, n_angles, n_h, n_v);

  numFunc++;
  fxk += real(0.5) * norm_pixels(tv2, n_angles, n_v, n_h);

  // BGS - does this duplicate at lot of the previous code suggesting a common
  // subroutine?
  int kk = -1;
  for (int k = 0; k <= k_max; k++) {
    add_output("Iteration ");
    add_output(k + 1);
    send_output();
    kk++;
    Lklist[kk] = bL;
    muklist[kk] = bmu;
    /* Calculate the gradient in yk */
    numGrad++;
    fyk = alpha * DTD(yk, Nablafyk, uijl, tau, Ddim, Dm, Dn, Dl, sz);

    scal_y(alpha, Nablafyk, nx, ny, nz);

    /*-----------------Forward projection--------------------------------*/
    init_data(tv2, n_angles, n_h, n_v);
    device->forward_project(tv2, yk, grid_offset, voxel_size, Dm, Dn, Dl);
    sum_axpy(real(-1.0), b, tv2, n_angles, n_h, n_v);

    numFunc++;
    fyk += real(0.5) * norm_pixels(tv2, n_angles, n_v, n_h);

    init_data(tv, nx, ny, nz);

    /*------------------Backward projection------------------------------*/
    init_data(temp, nx, ny, nz);
    device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
    /* For each voxel */
    sum_axpy(real(1.0), temp, Nablafyk, nx, ny,
	     nz);

    /* Update estimate of the strong convexity parameter as minimum of current
       value and the computed value between xk and yk */
    if (k != 0) {
      diff_xyz(xk, yk, tv, nx, ny, nz);

      bmu = std::max(std::min(real(2) * (fxk * real(one_eps)
					 - (fyk
					    + dot_prod(Nablafyk, tv, nx,ny,nz)))
			      / norm_voxels(tv, nx, ny, nz), bmu),
		     real(0.0));
    }

    /* Take the projected step from yk to xkp1 */
    t = - 1 / bL;
    sum_xbyz(yk, t, Nablafyk, xkp1, nx, ny, nz);

    P(xkp1, ctype, d, c, sz);

    /* Backtracking on Lipschitz parameter. */
    diff_xyz(xkp1, yk, tv, nx, ny, nz);

    hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl, sz);

    /*-----------------Forward projection--------------------------------*/
    init_data(tv2, n_angles, n_h, n_v);
    device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
    sum_axpy(real(-1.0), b, tv2, n_angles, n_h, n_v);

    gxkp1 = real(0.5) * norm_pixels(tv2, n_angles, n_v, n_h);

    numFunc++;
    fxkp1 = hxkp1 + gxkp1;

    while (fxkp1 / one_eps > fyk + dot_prod(Nablafyk, tv, nx, ny, nz)
	   + (bL / real(2)) * norm_voxels(tv, nx, ny, nz)) {
      numBack++;
      bL *= s_L;

      /* Take the projected step from yk to xkp1 */
      t = - 1 / bL;
      sum_xbyz(yk, t, Nablafyk, xkp1, nx, ny, nz);
      P(xkp1, ctype, d, c, sz);

      /* Backtracking on Lipschitz parameter. */
      diff_xyz(xkp1, yk, tv, nx, ny, nz);

      hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl, sz);

      /*-----------------Forward projection--------------------------------*/
      init_data(tv2, n_angles, n_h, n_v);
      device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
      sum_axpy(real(-1.0), b, tv2, n_angles, n_h, n_v);

      gxkp1 = real(0.5) * norm_pixels(tv2, n_angles, n_v, n_h);

      numFunc++;
      fxkp1 = hxkp1 + gxkp1;
    }

    fxkp1l[kk] = fxkp1;

    if (ghxl) {
      hxkp1l[kk] = hxkp1;
      gxkp1l[kk] = gxkp1;
    }

    /* store the iterate if requested */
    if (xl) {
      sl_int l = 0;
      for (int i = 0; i < int(sz[2]); i++) {
	for (int j = 0; j < int(sz[1]); j++) {
	  for (int k = 0; k < int(sz[0]); k++) {
	    xlist[kk * prodDims + l] = xkp1[k][j][i];
	    l++;
	  }
	}
      }
    }

    if(verbose)
      printf("k=%6d  f(x^k+1)=%e  ||G_L(x^k+1)||=%e  L_k=%.2e  mu_k=%.2e\n",
	     kk, fxkp1, nGt, bL, bmu);

    /* calculate the gradient in xkp1 */
    numGrad++;
    scal_y(alpha, Nablafxkp1, nx, ny, nz);

    /*------------------Backward projection------------------------------*/
    init_data(temp, nx, ny, nz);
    device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
    sum_axpy(real(1.0), temp, Nablafxkp1, nx, ny, nz);

    /* Check stopping criteria xkp1*/
    if (ctype == 1) {
      nGt = std::sqrt(norm_voxels(Nablafxkp1, nx, ny, nz));
      if (nGt <= epsb_rel * prodDims) {
	stop = true;
	/*overwrite xkp1 to return with*/
	t = - 1 / bL;
	sum_axpy(t, Nablafxkp1, xkp1, nx, ny, nz);
      }
    } else {
      t = - 1 / bL;
      sum_xbyz(xkp1, t, Nablafxkp1, tv, nx, ny, nz);
      P(tv, ctype, d, c, sz);

      scal_xby(xkp1, real(-1.0), tv, nx, ny, nz);

      nGt = bL * std::sqrt(norm_voxels(tv, nx, ny, nz));
      if (nGt <= epsb_rel * prodDims) {
	stop = true;
	/*overwrite xkp1 to return with*/
	sum_axpy(t, Nablafxkp1, xkp1, nx, ny, nz);
	P(xkp1, ctype, d, c, sz);
      }
    }

    if(stop or STOPMARK or kk == k_max) {
      goto cleanup;
    }

    /* Check stopping criteria yk*/
    if (ctype == 1) {
      if (std::sqrt(norm_voxels(Nablafyk, nx, ny, nz)) <= epsb_rel * prodDims)
	stop = true;
    } else {
      t = - 1 / bL;
      diff_xyz(yk, xkp1, tv, nx, ny, nz);
      if (bL * std::sqrt(norm_voxels(tv, nx, ny, nz)) <= epsb_rel * prodDims)
	stop = true;
    }

    if (stop or STOPMARK or kk == k_max)
      goto cleanup;

    /* Compute values for accelerated step size*/
    q = bmu / bL;

    if (k != 0)
      cumprod *= (real(1) - std::sqrt(q));
    /*tjeck if the convergence rate is fast enough*/
    if (bmu > 0) {
      if (nGt * nGt > cumprod * (real(4) * bL / bmu - bL / Lm1 + real(4) *gamma0
				 * bL / std::pow(bmu, 2)) * nGtm1 * nGtm1) {
	/*printf("not fast enough %d\n",kk);*/
	/*printf("%f %f %f %f\n",nGt,nGtm1*nGtm1,cumprod,
	       (4*bL/bmu-bL/Lm1+4*gamma0*bL/pow(bmu,2)));*/
	/*list of restart positions */
	rp.push_back(kk);
	bmu *= s_mu;
	numRest++;
	goto restart;
      }
    }
    /*printf("%f\n",q);DRAW;*/

    thetakp1 = (-(std::pow(thetak, 2) - q)
		+ std::sqrt(std::pow(std::pow(thetak, 2) - q, 2)
			    + real(4) * std::pow(thetak, 2))) / real(2.0);
    betak = (thetak * (real(1) - thetak)) / (std::pow(thetak, 2) + thetakp1);

    if (k == 0) {
      gamma0 = thetakp1 * (thetakp1 * bL - bmu) / (real(1) - thetakp1);
      /*printf("gamma0 %f\n",gamma0);*/
    }

    /* accelerated term*/
    /* yk = xkp1 + betak*(xkp1-xk) */
    diff_xyz(xkp1, xk, tv, nx, ny, nz);

    sum_xbyz(xkp1, betak, tv, yk, nx, ny, nz);

    /* Update the values for the next iteration*/
    thetak = thetakp1;
    copy(xkp1, xk, nx, ny, nz);

    fxk = fxkp1;

  }

 cleanup:
  //delete [] Nablafyk;
  //delete [] Nablafxkp1;

  //delete [] yk;
  //delete [] xk;

  //delete [] tv;
  //delete [] tv2;
  //delete [] uijl;
  //delete [] temp;

  kend = kk;
}

void P(voxel_data &y, const int ctype, const real_1dr &d,
       const real_1dr &c, const voxel_data::size_type sz[])
{
  if (ctype == 2) { /* c <= x <= d (elementwise) */
    sl_int mnl = 0;
    for (int i = 0; i < int(sz[2]); i++) {
      for (int j = 0; j < int(sz[1]); j++) {
	for (int k = 0; k < int(sz[0]); k++) {
	  if (y[k][j][i] < c[mnl])
	    y[k][j][i] = c[mnl];
	  else if (y[k][j][i] > d[mnl])
	    y[k][j][i] = d[mnl];
	}
	mnl++;
      }
    }
  } else if (ctype == 3) { /* c <= x <= d (elementwise) */
    for (int i = 0; i < int(sz[2]); i++) {
      for (int j = 0; j < int(sz[1]); j++) {
	for (int k = 0; k < int(sz[0]); k++) {
	  if (y[k][j][i] < c[0])
	    y[k][j][i] = c[0];
	  else if (y[k][j][i] > d[0])
	    y[k][j][i] = d[0];
	}
      }
    }
  }
}

/* Function used to calculate operations involving D and D^T*/
real DTD(voxel_data &x, voxel_data &Nablafx, std::vector<real> &uijl,
	 const real tau, const int Ddim, const int Dm, const int Dn,
	 const int Dl, const voxel_data::size_type sz[])
{
  real tv_tau_x = 0;
  real taud2 = tau / 2;
  real tau2 = 1 / (tau * 2);

  /* Clear the current gradient */
  init_data(Nablafx, sl_int(sz[0]), sl_int(sz[1]), sl_int(sz[2]));

  if (Ddim == 2) {
    for (sl_int u = 0; u < Dm; u++) {
      for (sl_int v = 0; v < Dn; v++) {
	sl_int i1 = (u + 1) % Dm;
	sl_int i2 = (v + 1) % Dn;

	uijl[0] = x[i1][v][0] - x[u][v][0];
	uijl[1] = x[u][i2][0] - x[u][v][0];

	real c1 = std::sqrt(uijl[0] * uijl[0] + uijl[1] * uijl[1]);
	real c2;
	if (tau < c1) {
	  c2 = c1;
	  tv_tau_x += c1 - taud2;
	} else {
	  c2 = tau;
	  tv_tau_x += c1 * c1 * tau2;
	}

	uijl[0] /= c2;
	uijl[1] /= c2;

	Nablafx[i1][v][0] += uijl[0];
	Nablafx[u][v][0] -= uijl[0];
	Nablafx[u][i2][0] += uijl[1];
	Nablafx[u][v][0] -= uijl[1];
      }
    }

  } else if (Ddim == 3) {
    for (sl_int u = 0; u < Dm; u++) {
      int u1 = (u + 1) % Dm;
      for (sl_int v = 0; v < Dn; v++) {
	int v1 = (v + 1) % Dn;
	//sl_int s1= ((u + 1) % Dm) + v * Dm;
	//sl_int s2 = u + ((v + 1) % Dn) * Dm;
	//sl_int s3 = u + v * Dm;
	//sl_int s4 = 0;
	for (sl_int w = 0; w < Dl; w++) {
	  int w1 = (w + 1) % Dl;
	  //sl_int i1 = s1 + s4;
	  //sl_int i2 = s2 + s4;
	  //sl_int i4 = s3 + s4;
	  //sl_int i3 = (i4 + mn) % mnl;

	  //s4 += mn;

	  uijl[0] = x[u1][v][w] - x[u][v][w];
	  uijl[1] = x[u][v1][w] - x[u][v][w];
	  uijl[2] = x[u][v][w1] - x[u][v][w];
	  //uijl[0] = (x.data())[i1] - (x.data())[i4];
	  //uijl[1] = (x.data())[i2] - (x.data())[i4];
	  //uijl[2] = (x.data())[i3] - (x.data())[i4];

	  real c1 = std::sqrt(uijl[0] * uijl[0] + uijl[1] * uijl[1]
			      + uijl[2] * uijl[2]);
	  real c2;
	  if (tau > c1) {
	    c2 = tau;
	    tv_tau_x += c1 * c1 * tau2;
	  } else {
	    c2 = c1;
	    tv_tau_x += c1 - taud2;
	  }

	  uijl[0] /= c2;
	  uijl[1] /= c2;
	  uijl[2] /= c2;

	  Nablafx[u1][v][w] += uijl[0];
	  Nablafx[u][v][w] -= uijl[0] + uijl[1] + uijl[2];
	  Nablafx[u][v1][w] += uijl[1];
	  Nablafx[u][v][w1] += uijl[2];
	  //(Nablafx.data())[i1] += uijl[0];
	  //(Nablafx.data())[i4] -= uijl[0] + uijl[1] + uijl[2];
	  //(Nablafx.data())[i2] += uijl[1];
	  //(Nablafx.data())[i3] += uijl[2];

	}
      }
    }
  } else
    report_error("Incorrect dim variable, only dim=2 or dim=3 supported.");

  return tv_tau_x;
}
