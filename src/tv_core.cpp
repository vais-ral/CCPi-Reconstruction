
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <algorithm>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // types
#include "blas.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "ui_calls.hpp"

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

static void P(voxel_data &y, const int ctype, const real d[], const real c[],
	      const voxel_data::size_type sz[]);
static real DTD(voxel_data &x, voxel_data &Nablafx, real uijl[],
		const real tau, const int Ddim, const int Dm, const int Dn,
		const int Dl, const voxel_data::size_type sz[]);

void CCPi::tv_regularization::tvreg_core(voxel_data &xkp1, real *fxkp1,
					 real *hxkp1, real *gxkp1, real *fxkp1l,
					 int *kend, const real voxel_size[],
					 const pixel_type *b, const real alpha,
					 real tau, real bL, real bmu,
					 real epsb_rel, int k_max,
					 const int Ddim, const int Dm,
					 const int Dn, const int Dl,
					 const sl_int prodDims, int ctype,
					 real *d, real *c, const bool ghxl,
					 const bool xl, real *hxkp1l,
					 real *gxkp1l, real *xlist,
					 const bool verbose, real *numGrad,
					 real *numBack, real *numFunc,
					 real *numRest, real *Lklist,
					 real *muklist, std::list<int> &rp,
					 const real grid_offset[],
					 const instrument *device)
{
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
  sl_int n_rays = device->get_data_size();

  INITBREAK

  const voxel_data::size_type *sz = xkp1.shape();
  //sl_int n_vox = sl_int(sz[0]) * sl_int(sz[1]) * sl_int(sz[2]);
  voxel_3d Nablafyk(boost::extents[sz[0]][sz[1]][sz[2]],
		    boost::fortran_storage_order());
  voxel_3d Nablafxkp1(boost::extents[sz[0]][sz[1]][sz[2]],
		      boost::fortran_storage_order());
  voxel_3d yk(boost::extents[sz[0]][sz[1]][sz[2]],
	      boost::fortran_storage_order());
  voxel_3d xk(boost::extents[sz[0]][sz[1]][sz[2]],
	      boost::fortran_storage_order());
  voxel_3d temp(boost::extents[sz[0]][sz[1]][sz[2]],
		boost::fortran_storage_order());

  /*temp vectors */
  voxel_3d tv(boost::extents[sz[0]][sz[1]][sz[2]],
	      boost::fortran_storage_order());
  pixel_type *tv2 = new pixel_type[n_rays];
  real *uijl = new real[Ddim];

  /* INITIALIZE */
  *numGrad = 0;
  *numBack = 0;
  *numFunc = 0;
  *numRest = 0;

 restart:
  real cumprod = 1.0;

  /* Project solution onto feasible space */
  P(xkp1, ctype, d, c, sz);

  real thetak = std::sqrt(bmu / bL);

  /* Include a backtracking to initialize everything */
  /* y_k = x_k+1 */
  dcopy(prodDims, xkp1, 1, yk, 1);

  /* Calculate the gradient in yk=xk+1 */
  (*numGrad)++;

  /* alpha*T_tau(y_k) (Nablafyk is returned as gradient) */
  real fyk = alpha * DTD(yk, Nablafyk, uijl, tau, Ddim, Dm, Dn, Dl, sz);

  /* For each voxel */
  for (int i = 0; i < int(sz[2]); i++)
    for (int j = 0; j < int(sz[1]); j++)
      for (int k = 0; k < int(sz[0]); k++)
	Nablafyk[k][j][i] *= alpha;
  /*-----------------Forward projection--------------------------------*/
  for (sl_int i = 0; i < n_rays; i++)
    tv2[i] = 0;
  /*tv2 = A*yk */
  device->forward_project(tv2, yk, grid_offset, voxel_size, Dm, Dn, Dl);
  /* tv2 = tv2 - b */
  for (sl_int i = 0; i < n_rays; i++)
    tv2[i] = tv2[i] - b[i];

  (*numFunc)++;
  /* fyk + 0.5*||A*y_k - b||^2 */
  fyk += real(0.5) * std::pow(dnrm2(n_rays, tv2, 1), 2);

  for (int i = 0; i < int(sz[2]); i++)
    for (int j = 0; j < int(sz[1]); j++)
      for (int k = 0; k < int(sz[0]); k++)
	tv[k][j][i] = 0.0;
  /*------------------Backward projection------------------------------*/
  for (int i = 0; i < int(sz[2]); i++)
    for (int j = 0; j < int(sz[1]); j++)
      for (int k = 0; k < int(sz[0]); k++)
	temp[k][j][i] = 0.0;
  device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
  /* For each voxel */
  for (int i = 0; i < int(sz[2]); i++)
    for (int j = 0; j < int(sz[1]); j++)
      for (int k = 0; k < int(sz[0]); k++)
	Nablafyk[k][j][i] += temp[k][j][i];

  /* Take the projected step from yk to xkp1 */
  /* bL is original setting of L_k */
  t = - 1 / bL;
  dcopy(prodDims, yk, 1, xkp1, 1);
  /* x_k+1 = y_k - t*Nablaf(y_k) */
  daxpy(prodDims, t, Nablafyk, 1, xkp1, 1);

  P(xkp1, ctype, d, c, sz);

  /* Backtracking on Lipschitz parameter. */
  for (int i = 0; i < int(sz[2]); i++)
    for (int j = 0; j < int(sz[1]); j++)
      for (int k = 0; k < int(sz[0]); k++)
	tv[k][j][i] = xkp1[k][j][i] - yk[k][j][i];

  /* alpha*T_tau(xkp1+1) (Nablafyk is returned as gradient) */
  *hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl, sz);

  /*-----------------Forward projection--------------------------------*/
  for (sl_int i = 0; i < n_rays; i++)
    tv2[i] = 0;
  /*tv2 = A*x_k+1 */
  device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
  /* tv2 = (A*x_k+1 - b) */
  for (sl_int i = 0; i < n_rays; i++)
    tv2[i] -= b[i];

  /* 0.5*||A*x_k+1 - b||^2 */
  *gxkp1 = real(0.5) * std::pow(dnrm2(n_rays, tv2, 1), 2);

  (*numFunc)++;
  // f(x_k+1) = h(x_k+1) + g(x_k+1) = alpha*T_tau(x_k+1) + 0.5*||A*x_k+1 - b||^2
  *fxkp1 = *hxkp1 + *gxkp1;

  while (*fxkp1 / (1+1e-14) > fyk + ddot(prodDims, Nablafyk, 1, tv, 1)
	 + (bL / 2) * std::pow(dnrm2(prodDims, tv, 1), 2)) {
    (*numBack)++;
    bL *= s_L;

    /* Take the projected step from yk to xkp1 */
    /* t = - L^-1 */
    t = - 1 / bL;
    /* x_k+1 = y_k */
    dcopy(prodDims, yk, 1, xkp1, 1);
    /* x_k+1 = x_k+1 - t*Nablafyk */
    daxpy(prodDims, t, Nablafyk, 1, xkp1, 1);
    /* project to within bounds (x_k+1) */
    P(xkp1, ctype, d, c, sz);

    /* Backtracking on Lipschitz parameter. */
    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  tv[k][j][i] = xkp1[k][j][i] - yk[k][j][i];

    /* alpha*T_tau(x_k+1) */
    *hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl, sz);

    /*-----------------Forward projection--------------------------------*/
    for (sl_int i = 0; i < n_rays; i++)
      tv2[i] = 0;
    /*tv2 = A*x_k+1 */
    device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
    /* tv2 = (A*x_k+1 - b) */
    for (sl_int i = 0; i < n_rays; i++)
      tv2[i] -= b[i];

    /* 0.5*||A*x_k+1 - b||^2 */
    *gxkp1 = real(0.5) * std::pow(dnrm2(n_rays, tv2, 1), 2);

    (*numFunc)++;
    /* f(x_k+1) = h(x_k+1) + g(x_k+1)
       = alpha*T_tau(x_k+1) + 0.5*||A*x_k+1 - b||^2 */
    *fxkp1 = *hxkp1 + *gxkp1;
  }

  /* Calculate initial gradient map */
  if (ctype == 1) { // No constraints - just real numbers
    nGt = dnrm2(prodDims, Nablafxkp1, 1);
  } else { /* Positivity constraint */
    t = - 1 / bL;
    dcopy(prodDims, xkp1, 1, tv, 1);
    daxpy(prodDims, t, Nablafxkp1, 1, tv, 1);
    P(tv, ctype, d, c, sz);

    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  tv[k][j][i] = xkp1[k][j][i] - tv[k][j][i];

    nGt = bL * dnrm2(prodDims, tv, 1);
  }

  /* save the initial parameter */
  Lm1 = bL;
  nGtm1 = nGt;
  dcopy(prodDims, xkp1, 1, xk, 1);
  dcopy(prodDims, xkp1, 1, yk, 1);

  /* LOOP */
  bool stop = false; /*Flag for when to break the for-loop*/

  /*Calculate fxk */
  real fxk = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,Ddim,Dm,Dn,Dl, sz);

  /*-----------------Forward projection--------------------------------*/
  for (sl_int i = 0; i < n_rays; i++)
    tv2[i] = 0;
  device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
  for (sl_int i = 0; i < n_rays; i++)
    tv2[i] -= b[i];

  (*numFunc)++;
  fxk += real(0.5) * std::pow(dnrm2(n_rays, tv2, 1), 2);

  // BGS - does this duplicate at lot of the previous code suggesting a common
  // subroutine?
  int kk = -1;
  for (int k = 0; k <= k_max; k++) {
    kk++;
    Lklist[kk] = bL;
    muklist[kk] = bmu;
    /* Calculate the gradient in yk */
    (*numGrad)++;
    fyk = alpha * DTD(yk, Nablafyk, uijl, tau, Ddim, Dm, Dn, Dl, sz);

    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  Nablafyk[k][j][i] *= alpha;

    /*-----------------Forward projection--------------------------------*/
    for (sl_int i = 0; i < n_rays; i++)
      tv2[i] = 0;
    device->forward_project(tv2, yk, grid_offset, voxel_size, Dm, Dn, Dl);
    for (sl_int i = 0; i < n_rays; i++)
      tv2[i] -= b[i];

    (*numFunc)++;
    fyk += real(0.5) * std::pow(dnrm2(n_rays, tv2, 1), 2);

    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  tv[k][j][i] = 0.0;

    /*------------------Backward projection------------------------------*/
    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  temp[k][j][i] = 0.0;
    device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
    /* For each voxel */
    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  Nablafyk[k][j][i] += temp[k][j][i];

    /* Update estimate of the strong convexity parameter as minimum of current
       value and the computed value between xk and yk */
    if (k != 0) {
      for (int i = 0; i < int(sz[2]); i++)
	for (int j = 0; j < int(sz[1]); j++)
	  for (int k = 0; k < int(sz[0]); k++)
	    tv[k][j][i] = xk[k][j][i] - yk[k][j][i];

      bmu = std::max(std::min(real(2) * (fxk * (real(1) + real(1e-14))
					 - (fyk + ddot(prodDims, Nablafyk,
						       1, tv, 1)))
			      / std::pow(dnrm2(prodDims, tv, 1), 2), bmu),
		     real(0.0));
    }

    /* Take the projected step from yk to xkp1 */
    t = - 1 / bL;
    dcopy(prodDims, yk, 1, xkp1, 1);
    daxpy(prodDims, t, Nablafyk, 1, xkp1, 1);

    P(xkp1, ctype, d, c, sz);

    /* Backtracking on Lipschitz parameter. */
    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  tv[k][j][i] = xkp1[k][j][i] - yk[k][j][i];

    *hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl, sz);

    /*-----------------Forward projection--------------------------------*/
    for (sl_int i = 0; i < n_rays; i++)
      tv2[i] = 0;
    device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
    for (sl_int i = 0; i < n_rays; i++)
      tv2[i] -= b[i];

    *gxkp1 = real(0.5) * std::pow(dnrm2(n_rays, tv2, 1), 2);

    (*numFunc)++;
    *fxkp1 = *hxkp1 + *gxkp1;

    while (*fxkp1 / (real(1) + real(1e-14)) > fyk
	   + ddot(prodDims, Nablafyk, 1, tv, 1)
	   + (bL / real(2)) * std::pow(dnrm2(prodDims, tv, 1), 2)) {
      (*numBack)++;
      bL *= s_L;

      /* Take the projected step from yk to xkp1 */
      t = - 1 / bL;
      dcopy(prodDims, yk, 1, xkp1, 1);
      daxpy(prodDims, t, Nablafyk, 1, xkp1, 1);
      P(xkp1, ctype, d, c, sz);

      /* Backtracking on Lipschitz parameter. */
      for (int i = 0; i < int(sz[2]); i++)
	for (int j = 0; j < int(sz[1]); j++)
	  for (int k = 0; k < int(sz[0]); k++)
	    tv[k][j][i] = xkp1[k][j][i] - yk[k][j][i];

      *hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl, sz);

      /*-----------------Forward projection--------------------------------*/
      for (sl_int i = 0; i < n_rays; i++)
	tv2[i] = 0;
      device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
      for (sl_int i = 0; i < n_rays; i++)
	tv2[i] -= b[i];

      *gxkp1 = real(0.5) * std::pow(dnrm2(n_rays, tv2, 1), 2);

      (*numFunc)++;
      *fxkp1 = *hxkp1 + *gxkp1;
    }

    fxkp1l[kk] = *fxkp1;

    if (ghxl) {
      hxkp1l[kk] = *hxkp1;
      gxkp1l[kk] = *gxkp1;
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
	     kk, *fxkp1, nGt, bL, bmu);

    /* calculate the gradient in xkp1 */
    (*numGrad)++;
    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  Nablafxkp1[k][j][i] *= alpha;

    /*------------------Backward projection------------------------------*/
    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  temp[k][j][i] = 0;
    device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  Nablafxkp1[k][j][i] += temp[k][j][i];

    /* Check stopping criteria xkp1*/
    if (ctype == 1) {
      nGt = dnrm2(prodDims, Nablafxkp1, 1);
      if (nGt <= epsb_rel * prodDims) {
	stop = true;
	/*overwrite xkp1 to return with*/
	t = - 1 / bL;
	daxpy(prodDims, t, Nablafxkp1, 1, xkp1, 1);
      }
    } else {
      t = - 1 / bL;
      dcopy(prodDims, xkp1, 1, tv, 1);
      daxpy(prodDims, t, Nablafxkp1, 1, tv, 1);
      P(tv, ctype, d, c, sz);

      for (int i = 0; i < int(sz[2]); i++)
	for (int j = 0; j < int(sz[1]); j++)
	  for (int k = 0; k < int(sz[0]); k++)
	    tv[k][j][i] = xkp1[k][j][i] - tv[k][j][i];

      nGt = bL * dnrm2(prodDims, tv, 1);
      if (nGt <= epsb_rel * prodDims) {
	stop = true;
	/*overwrite xkp1 to return with*/
	daxpy(prodDims, t, Nablafxkp1, 1, xkp1, 1);
	P(xkp1, ctype, d, c, sz);
      }
    }

    if(stop or STOPMARK or kk == k_max) {
      goto cleanup;
    }

    /* Check stopping criteria yk*/
    if (ctype == 1) {
      if (dnrm2(prodDims, Nablafyk, 1) <= epsb_rel * prodDims)
	stop = true;
    } else {
      t = - 1 / bL;
      for (int i = 0; i < int(sz[2]); i++)
	for (int j = 0; j < int(sz[1]); j++)
	  for (int k = 0; k < int(sz[0]); k++)
	    tv[k][j][i] = yk[k][j][i] - xkp1[k][j][i];
      if (bL * dnrm2(prodDims, tv, 1) <= epsb_rel * prodDims)
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
	(*numRest)++;
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
    for (int i = 0; i < int(sz[2]); i++)
      for (int j = 0; j < int(sz[1]); j++)
	for (int k = 0; k < int(sz[0]); k++)
	  tv[k][j][i] = xkp1[k][j][i] - xk[k][j][i];

    dcopy(prodDims, xkp1, 1, yk, 1);
    daxpy(prodDims, betak, tv, 1, yk, 1);

    /* Update the values for the next iteration*/
    thetak = thetakp1;
    dcopy(prodDims, xkp1, 1, xk, 1);

    fxk = *fxkp1;

  }

 cleanup:
  //delete [] Nablafyk;
  //delete [] Nablafxkp1;

  //delete [] yk;
  //delete [] xk;

  //delete [] tv;
  delete [] tv2;
  delete [] uijl;
  //delete [] temp;

  *kend = kk;
}

void P(voxel_data &y, const int ctype, const real d[], const real c[],
       const voxel_data::size_type sz[])
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
	  if (y[k][j][i] < *c)
	    y[k][j][i] = *c;
	  else if (y[k][j][i] > *d)
	    y[k][j][i] = *d;
	}
      }
    }
  }
}

/* Function used to calculate operations involving D and D^T*/
real DTD(voxel_data &x, voxel_data &Nablafx, real uijl[], const real tau,
	 const int Ddim, const int Dm, const int Dn, const int Dl,
	 const voxel_data::size_type sz[])
{
  real tv_tau_x=0;
  real taud2 = tau / 2;
  real tau2 = 1 / (tau * 2);

  /* Clear the current gradient */
  for (int i = 0; i < int(sz[2]); i++)
    for (int j = 0; j < int(sz[1]); j++)
      for (int k = 0; k < int(sz[0]); k++)
	Nablafx[k][j][i] = 0.0;

  if (Ddim == 2) {
    for (sl_int u = 0; u <= Dm - 1; u++) {
      for (sl_int v = 0; v <= Dn - 1; v++) {
	sl_int i1 = (u + 1) %Dm + v * Dm;
	sl_int i2 = u + ((v + 1) % Dn) * Dm;
	sl_int i3= u + v * Dm;

	uijl[0] = (x.data())[i1] - (x.data())[i3];
	uijl[1] = (x.data())[i2] - (x.data())[i3];

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

	(Nablafx.data())[i1] += uijl[0];
	(Nablafx.data())[i3] -= uijl[0];
	(Nablafx.data())[i2] += uijl[1];
	(Nablafx.data())[i3] -= uijl[1];
      }
    }

  } else if (Ddim == 3) {
    sl_int mn = Dm * Dn;
    sl_int mnl = mn * Dl;

    for (sl_int u = 0; u <= Dm - 1; u++) {
      for (sl_int v = 0; v <= Dn - 1; v++) {
	/*s1= ((u+1)%Dm) + v*Dm;*/
	/*s2 = u + ((v+1)%Dn)*Dm;*/
	sl_int s1= ((u + 1) % Dm) + v * Dm;
	sl_int s2 = u + ((v + 1) % Dn) * Dm;
	sl_int s3 = u + v * Dm;
	sl_int s4 = 0;
	for (sl_int w = 0; w <= Dl - 1; w++) {
	  sl_int i1 = s1 + s4;
	  sl_int i2 = s2 + s4;
	  sl_int i4 = s3 + s4;
	  sl_int i3 = (i4 + mn) % mnl;

	  s4 += mn;

	  uijl[0] = (x.data())[i1] - (x.data())[i4];
	  uijl[1] = (x.data())[i2] - (x.data())[i4];
	  uijl[2] = (x.data())[i3] - (x.data())[i4];

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

	  (Nablafx.data())[i1] += uijl[0];
	  (Nablafx.data())[i4] -= uijl[0] + uijl[1] + uijl[2];
	  (Nablafx.data())[i2] += uijl[1];
	  (Nablafx.data())[i3] += uijl[2];

	}
      }
    }
  } else
    report_error("Incorrect dim variable, only dim=2 or dim=3 supported.");

  return tv_tau_x;
}
