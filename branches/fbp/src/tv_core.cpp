
#include <cmath>
#include <cstdlib>
#include <cstdio>
#ifdef MATLAB_MEX_FILE
#  include "mex_types.hpp"
#else
#  include "base_types.hpp"
#endif // types
#include "blas.hpp"
#include "fbp.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"

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

static void P(voxel_type *y, const int ctype, const real d[],
	      const real c[], const long mnl);
static real DTD(voxel_type x[], voxel_type Nablafx[], real uijl[],
		const real tau, const int Ddim, const int Dm, const int Dn,
		const int Dl, const long prodDims);

void CCPi::tvreg_core(voxel_type *xkp1, real *fxkp1, real *hxkp1, real *gxkp1,
		      real *fxkp1l, int *kend, const real voxel_size[],
		      const real *b, const real alpha, real tau, real bL,
		      real bmu, real epsb_rel, int k_max, const int Ddim,
		      const int Dm, const int Dn, const int Dl,
		      const long prodDims, int ctype, real *d, real *c,
		      const bool ghxl, const bool xl, real *hxkp1l,
		      real *gxkp1l, real *xlist, const bool verbose,
		      real *numGrad, real *numBack, real *numFunc,
		      real *numRest, real *Lklist, real *muklist,
		      std::list<int> &rp, const real grid_offset[],
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
  long n_rays = device->get_data_size();

  INITBREAK

  voxel_type *Nablafyk = new voxel_type[prodDims];
  voxel_type *Nablafxkp1 = new voxel_type[prodDims];
  voxel_type *yk = new voxel_type[prodDims];
  voxel_type *xk = new voxel_type[prodDims];
  voxel_type *temp = new voxel_type[prodDims];

  /*temp vectors */
  voxel_type *tv = new voxel_type[prodDims];
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
  P(xkp1, ctype, d, c, prodDims);

  real thetak = std::sqrt(bmu / bL);

  /* Include a backtracking to initialize everything */
  /* y_k = x_k+1 */
  dcopy(prodDims, xkp1, 1, yk, 1);

  /* Calculate the gradient in yk=xk+1 */
  (*numGrad)++;

  /* alpha*T_tau(y_k) (Nablafyk is returned as gradient) */
  real fyk = alpha * DTD(yk, Nablafyk, uijl, tau, Ddim, Dm, Dn, Dl, prodDims);

  /* For each voxel */
  for (long i = 0; i < prodDims; i++)
    Nablafyk[i] *= alpha;
  /*-----------------Forward projection--------------------------------*/
  for (long i = 0; i < n_rays; i++)
    tv2[i] = 0;
  /*tv2 = A*yk */
  device->forward_project(tv2, yk, grid_offset, voxel_size, Dm, Dn, Dl);
  /* tv2 = tv2 - b */
  for (long i = 0; i < n_rays; i++)
    tv2[i] = tv2[i] - b[i];

  (*numFunc)++;
  /* fyk + 0.5*||A*y_k - b||^2 */
  fyk += 0.5 * pow(dnrm2(n_rays, tv2, 1), 2);

  for (long i = 0; i < prodDims; i++)
    tv[i] = 0;
  /*------------------Backward projection------------------------------*/
  for (long i = 0; i < prodDims; i++)
    temp[i] = 0;
  device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
  /* For each voxel */
  for (long i = 0; i < prodDims; i++)
    Nablafyk[i] += temp[i];

  /* Take the projected step from yk to xkp1 */
  /* bL is original setting of L_k */
  t = - 1 / bL;
  dcopy(prodDims, yk, 1, xkp1, 1);
  /* x_k+1 = y_k - t*Nablaf(y_k) */
  daxpy(prodDims, t, Nablafyk, 1, xkp1, 1);

  P(xkp1, ctype, d, c, prodDims);

  /* Backtracking on Lipschitz parameter. */
  for (long i = 0; i < prodDims; i++)
    tv[i] = xkp1[i]-yk[i];

  /* alpha*T_tau(xkp1+1) (Nablafyk is returned as gradient) */
  *hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl,
		       prodDims);

  /*-----------------Forward projection--------------------------------*/
  for (long i = 0; i < n_rays; i++)
    tv2[i] = 0;
  /*tv2 = A*x_k+1 */
  device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
  /* tv2 = (A*x_k+1 - b) */
  for (long i = 0; i < n_rays; i++)
    tv2[i] -= b[i];

  /* 0.5*||A*x_k+1 - b||^2 */
  *gxkp1 = 0.5 * pow(dnrm2(n_rays, tv2, 1), 2);

  (*numFunc)++;
  // f(x_k+1) = h(x_k+1) + g(x_k+1) = alpha*T_tau(x_k+1) + 0.5*||A*x_k+1 - b||^2
  *fxkp1 = *hxkp1 + *gxkp1;

  while (*fxkp1 / (1+1e-14) > fyk + ddot(prodDims, Nablafyk, 1, tv, 1)
	 + (bL / 2) * pow(dnrm2(prodDims, tv, 1), 2)) {
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
    P(xkp1, ctype, d, c, prodDims);

    /* Backtracking on Lipschitz parameter. */
    for (long i = 0; i < prodDims; i++)
      tv[i] = xkp1[i] - yk[i];

    /* alpha*T_tau(x_k+1) */
    *hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl,
			 prodDims);

    /*-----------------Forward projection--------------------------------*/
    for (long i = 0; i < n_rays; i++)
      tv2[i] = 0;
    /*tv2 = A*x_k+1 */
    device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
    /* tv2 = (A*x_k+1 - b) */
    for (long i = 0; i < n_rays; i++)
      tv2[i] -= b[i];

    /* 0.5*||A*x_k+1 - b||^2 */
    *gxkp1 = 0.5 * pow(dnrm2(n_rays, tv2, 1), 2);

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
    P(tv, ctype, d, c, prodDims);

    for (long i = 0; i < prodDims; i++)
      tv[i] = xkp1[i] - tv[i];

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
  real fxk = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,Ddim,Dm,Dn,Dl, prodDims);

  /*-----------------Forward projection--------------------------------*/
  for (long i = 0; i < n_rays; i++)
    tv2[i] = 0;
  device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
  for (long i = 0; i < n_rays; i++)
    tv2[i] -= b[i];

  (*numFunc)++;
  fxk += 0.5 * pow(dnrm2(n_rays, tv2, 1), 2);

  // BGS - does this duplicate at lot of the previous code suggesting a common
  // subroutine?
  int kk = -1;
  for (int k = 0; k <= k_max; k++) {
    kk++;
    Lklist[kk] = bL;
    muklist[kk] = bmu;
    /* Calculate the gradient in yk */
    (*numGrad)++;
    fyk = alpha * DTD(yk, Nablafyk, uijl, tau, Ddim, Dm, Dn, Dl, prodDims);

    for (long i = 0; i < prodDims; i++)
      Nablafyk[i] *= alpha;

    /*-----------------Forward projection--------------------------------*/
    for (long i = 0; i < n_rays; i++)
      tv2[i] = 0;
    device->forward_project(tv2, yk, grid_offset, voxel_size, Dm, Dn, Dl);
    for (long i = 0; i < n_rays; i++)
      tv2[i] -= b[i];

    (*numFunc)++;
    fyk += 0.5 * pow(dnrm2(n_rays, tv2, 1), 2);

    for (long i = 0; i < prodDims; i++)
      tv[i] = 0;

    /*------------------Backward projection------------------------------*/
    for (long i = 0; i < prodDims; i++)
      temp[i] = 0;
    device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
    /* For each voxel */
    for (long i = 0; i < prodDims; i++)
      Nablafyk[i] += temp[i];

    /* Update estimate of the strong convexity parameter as minimum of current
       value and the computed value between xk and yk */
    if (k != 0) {
      for (long i = 0; i < prodDims; i++)
	tv[i] = xk[i] - yk[i];

      bmu = std::max(std::min(2 * (fxk * (1 + 1e-14)
				   - (fyk + ddot(prodDims, Nablafyk, 1, tv, 1)))
			      / pow(dnrm2(prodDims, tv, 1), 2), bmu), 0.0);
    }

    /* Take the projected step from yk to xkp1 */
    t = - 1 / bL;
    dcopy(prodDims, yk, 1, xkp1, 1);
    daxpy(prodDims, t, Nablafyk, 1, xkp1, 1);

    P(xkp1, ctype, d, c, prodDims);

    /* Backtracking on Lipschitz parameter. */
    for (long i = 0; i < prodDims; i++)
      tv[i] = xkp1[i] - yk[i];

    *hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl,
			 prodDims);

    /*-----------------Forward projection--------------------------------*/
    for (long i = 0; i < n_rays; i++)
      tv2[i] = 0;
    device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
    for (long i = 0; i < n_rays; i++)
      tv2[i] -= b[i];

    *gxkp1 = 0.5 * pow(dnrm2(n_rays, tv2, 1), 2);

    (*numFunc)++;
    *fxkp1 = *hxkp1 + *gxkp1;

    while (*fxkp1 / (1 + 1e-14) > fyk + ddot(prodDims, Nablafyk, 1, tv, 1)
	   + (bL / 2) * pow(dnrm2(prodDims, tv, 1), 2)) {
      (*numBack)++;
      bL *= s_L;

      /* Take the projected step from yk to xkp1 */
      t = - 1 / bL;
      dcopy(prodDims, yk, 1, xkp1, 1);
      daxpy(prodDims, t, Nablafyk, 1, xkp1, 1);
      P(xkp1, ctype, d, c, prodDims);

      /* Backtracking on Lipschitz parameter. */
      for (long i = 0; i < prodDims; i++)
	tv[i] = xkp1[i] - yk[i];

      *hxkp1 = alpha * DTD(xkp1, Nablafxkp1, uijl, tau, Ddim, Dm, Dn, Dl,
			   prodDims);

      /*-----------------Forward projection--------------------------------*/
      for (long i = 0; i < n_rays; i++)
	tv2[i] = 0;
      device->forward_project(tv2, xkp1, grid_offset, voxel_size, Dm, Dn, Dl);
      for (long i = 0; i < n_rays; i++)
	tv2[i] -= b[i];

      *gxkp1 = 0.5 * pow(dnrm2(n_rays, tv2, 1), 2);

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
      for (long i = 0; i < prodDims; i++) {
	xlist[kk * prodDims + i] = xkp1[i];
      }
    }

    if(verbose)
      printf("k=%6d  f(x^k+1)=%e  ||G_L(x^k+1)||=%e  L_k=%.2e  mu_k=%.2e\n",
	     kk, *fxkp1, nGt, bL, bmu);

    /* calculate the gradient in xkp1 */
    (*numGrad)++;
    for (long i = 0; i < prodDims; i++)
      Nablafxkp1[i] *= alpha;

    /*------------------Backward projection------------------------------*/
    for (long i = 0; i < prodDims; i++)
      temp[i] = 0;
    device->backward_project(tv2, temp, grid_offset, voxel_size, Dm, Dn, Dl);
    for (long i = 0; i < prodDims; i++)
      Nablafxkp1[i] += temp[i];

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
      P(tv, ctype, d, c, prodDims);

      for (long i = 0; i < prodDims; i++)
	tv[i] = xkp1[i] - tv[i];

      nGt = bL * dnrm2(prodDims, tv, 1);
      if (nGt <= epsb_rel * prodDims) {
	stop = true;
	/*overwrite xkp1 to return with*/
	daxpy(prodDims, t, Nablafxkp1, 1, xkp1, 1);
	P(xkp1, ctype, d, c, prodDims);
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
      for (long i = 0; i < prodDims; i++)
	tv[i] = yk[i] - xkp1[i];
      if (bL * dnrm2(prodDims, tv, 1) <= epsb_rel * prodDims)
	stop = true;
    }

    if (stop or STOPMARK or kk == k_max)
      goto cleanup;

    /* Compute values for accelerated step size*/
    q = bmu / bL;

    if (k != 0)
      cumprod *= (1 - std::sqrt(q));
    /*tjeck if the convergence rate is fast enough*/
    if (bmu > 0) {
      if (nGt * nGt > cumprod * (4 * bL / bmu - bL / Lm1 + 4 *gamma0
				 * bL / pow(bmu, 2)) * nGtm1 * nGtm1) {
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

    thetakp1 = (-(pow(thetak, 2) - q) + std::sqrt(pow(pow(thetak, 2) - q, 2)
						  + 4 * pow(thetak, 2))) / 2.0;
    betak = (thetak * (1 - thetak)) / (pow(thetak, 2) + thetakp1);

    if (k == 0) {
      gamma0 = thetakp1 * (thetakp1 * bL - bmu) / (1 - thetakp1);
      /*printf("gamma0 %f\n",gamma0);*/
    }

    /* accelerated term*/
    /* yk = xkp1 + betak*(xkp1-xk) */
    for (long i = 0; i < prodDims; i++)
      tv[i] = xkp1[i] - xk[i];

    dcopy(prodDims, xkp1, 1, yk, 1);
    daxpy(prodDims, betak, tv, 1, yk, 1);

    /* Update the values for the next iteration*/
    thetak = thetakp1;
    dcopy(prodDims, xkp1, 1, xk, 1);

    fxk = *fxkp1;

  }

 cleanup:
  delete [] Nablafyk;
  delete [] Nablafxkp1;

  delete [] yk;
  delete [] xk;

  delete [] tv;
  delete [] tv2;
  delete [] uijl;
  delete [] temp;

  *kend = kk;
}

/* Functions used for the inverse problems using Nesterov or BB method
  Project onto the feasible (convex) set.
   c==1: Unconstrained.
   c==2: Lower and upper bounds (elementwise) on x. Inplace.
   c==3: as 2 but single value bounds applied whole space. Inplace. */
void P(voxel_type *y, const int ctype, const real d[], const real c[],
       const long mnl)
{
  if (ctype == 2) { /* c <= x <= d (elementwise) */
    for (long i = 0; i < mnl; i++) {
      if (y[i] < c[i])
	y[i] = c[i];
      else if (y[i] > d[i])
	y[i] = d[i];
    }
  } else if (ctype == 3) { /* c <= x <= d (elementwise) */
    for (long i = 0; i < mnl; i++) {
      if (y[i] < *c)
	y[i] = *c;
      else if (y[i] > *d)
	y[i] = *d;
    }
  }
}

/* Function used to calculate operations involving D and D^T*/
real DTD(voxel_type x[], voxel_type Nablafx[], real uijl[], const real tau,
	 const int Ddim, const int Dm, const int Dn, const int Dl,
	 const long prodDims)
{
  real tv_tau_x=0;
  real taud2 = tau / 2;
  real tau2 = 1 / (tau * 2);

  /* Clear the current gradient */
  for (long i = 0; i < prodDims; i++)
    Nablafx[i] = 0.0;

  if (Ddim == 2) {
    for (long u = 0; u <= Dm - 1; u++) {
      for (long v = 0; v <= Dn - 1; v++) {
	long i1 = (u + 1) %Dm + v * Dm;
	long i2 = u + ((v + 1) % Dn) * Dm;
	long i3= u + v * Dm;

	uijl[0] = x[i1] - x[i3];
	uijl[1] = x[i2] - x[i3];

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

	Nablafx[i1] += uijl[0];
	Nablafx[i3] -= uijl[0];
	Nablafx[i2] += uijl[1];
	Nablafx[i3] -= uijl[1];
      }
    }

  } else if (Ddim == 3) {
    long mn = Dm * Dn;
    long mnl = mn * Dl;

    for (long u = 0; u <= Dm - 1; u++) {
      for (long v = 0; v <= Dn - 1; v++) {
	/*s1= ((u+1)%Dm) + v*Dm;*/
	/*s2 = u + ((v+1)%Dn)*Dm;*/
	long s1= ((u + 1) % Dm) + v * Dm;
	long s2 = u + ((v + 1) % Dn) * Dm;
	long s3 = u + v * Dm;
	long s4 = 0;
	for (long w = 0; w <= Dl - 1; w++) {
	  long i1 = s1 + s4;
	  long i2 = s2 + s4;
	  long i4 = s3 + s4;
	  long i3 = (i4 + mn) % mnl;

	  s4 += mn;

	  uijl[0] = x[i1] - x[i4];
	  uijl[1] = x[i2] - x[i4];
	  uijl[2] = x[i3] - x[i4];

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

	  Nablafx[i1] += uijl[0];
	  Nablafx[i4] -= uijl[0] + uijl[1] + uijl[2];
	  Nablafx[i2] += uijl[1];
	  Nablafx[i3] += uijl[2];

	}
      }
    }
  } else
    printf("Incorrect dim variable, only dim=2 or dim=3 supported.\n");

  return tv_tau_x;
}
