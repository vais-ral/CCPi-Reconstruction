#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "tools.h"
#include "jacobs_rays.h"

/* Settings which makes the user do a CTRL-C break out of the loop*/
#if defined(LIBUT) && (defined(_WIN32) || defined(__WIN32__) )

#define STOPMARK utIsInterruptPending()
#define INITBREAK ;
bool utIsInterruptPending(void);

#else

#include <signal.h>
#define INITBREAK   sigint_cont = 1;  (void) signal(SIGINT , ex_sigint);
#define STOPMARK sigint_cont==0
int sigint_cont = 1;
void ex_sigint(int sig) {
	sigint_cont = 0;
}
#endif

extern void forwardProjection(double *source_x, double *source_y,
			      double *source_z, double *det_x, double *det_y,
			      double *det_z, double *ray_data, double *vol_data,
			      double *angles, struct jacobs_options *options,
			      int n_angles, long n_rays_y, long n_rays_z);
extern void backwardProjection(double *source_x, double *source_y,
			       double *source_z, double *det_x, double *det_y,
			       double *det_z, double *ray_data, double *vol_data,
			       double *angles, struct jacobs_options *options,
			       int n_angles, long n_rays_y, long n_rays_z,
			       int size_z, double grid_offset[],
			       double voxel_size[]);

void tvreg_upn_core(double *xkp1,double *fxkp1,double *hxkp1,double *gxkp1,double *fxkp1l,double *kend, double *voxel_size, double *dims, double *b,double alpha,double tau,
                    double bL,double bmu,double epsb_rel,int k_max,Dtype D,int ctype,double *d,double *c, int ghxl, int xl,double *hxkp1l,double *gxkp1l,double *xlist,
                    int verbose,double *numGrad,double* numBack,double *numFunc,double *numRest,double *Lklist,double *muklist,listelement *p_rp,
                    double *source_x, double *source_y, double *source_z, double *det_x, double *det_y, double *det_z, double *angles, double *grid_offset,
                    int n_rays_y, int n_rays_z, int n_angles){

  printf("IN TVREG_UPN_CORE \n");  
    
  double *yk,*xk,*Nablafyk,*Nablafxkp1,*uijl,*tv,*tv2, *temp;

  double fxk, fyk, q, thetak, thetakp1, betak, nGt, t, s_L=1.3, s_mu=0.7, Lm1, nGtm1, gamma0, cumprod;
  int i, j, k=0, prodDims, stop=0, one=1, kk=-1;
  

  struct jacobs_options options;
  int n_rays;

  /* set parameters for projections */
  
  options.im_size_default = 0;
  options.im_size_x = dims[0];
  options.im_size_y = dims[1];
  options.im_size_z = dims[2];

  options.b_default = 0;
  options.b_x = grid_offset[0];
  options.b_y = grid_offset[1];
  options.b_z = grid_offset[2];

  options.d_default = 0;
  options.d_x = voxel_size[0];
  options.d_y = voxel_size[1]; 
  options.d_z = voxel_size[2];

  n_rays = n_rays_y*n_rays_z*n_angles;

  INITBREAK

	prodDims=D.prodDims;
  
	Nablafyk = malloc(prodDims*sizeof(double));
	Nablafxkp1 = malloc(prodDims*sizeof(double));
	yk = malloc(prodDims*sizeof(double));
	xk = malloc(prodDims*sizeof(double));
	temp = malloc(prodDims*sizeof(double));
    
	/*temp vectors */
	tv = malloc(prodDims*sizeof(double));
    	tv2 = malloc(n_rays*sizeof(double));
   	uijl = malloc(D.dim*sizeof(double));
  
	
  /* INITIALIZE */
	numGrad[0]=0;numBack[0]=0;numFunc[0]=0;numRest[0]=0;
   
 restart:
	cumprod=1.0;

  	/* Project solution onto feasible space */
	P(xkp1,ctype,d,c,prodDims);

	thetak=sqrt(bmu/bL);

	/* Include a backtracking to initialize everything */

	dcopy_(&prodDims,xkp1,&one,yk,&one); /* y_k = x_k+1 */

	/* Calculate the gradient in yk=xk+1 */							 
	numGrad[0]+=1;

	fyk = alpha*DTD(yk,Nablafyk,uijl,tau,D); /* alpha*T_tau(y_k) (Nablafyk is returned as gradient) */

	for(i=0;i<prodDims;i++){	Nablafyk[i]=alpha*Nablafyk[i]; } /* For each voxel */
	
	/*-----------------Forward projection--------------------------------*/
	for(i=0;i<n_rays;i++){	tv2[i]=0; }
    	forwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, yk, angles, &options, n_angles, n_rays_y, n_rays_z); /*tv2 = A*yk */
	for(i=0;i<n_rays;i++){	tv2[i] = tv2[i]-b[i]; } /* tv2 = tv2 - b */

        numFunc[0]+=1;
	fyk = fyk + 0.5*pow( dnrm2_(&n_rays,tv2,&one),2);  /* fyk + 0.5*||A*y_k - b||^2 */

	for(i=0;i<prodDims;i++){	tv[i]=0; }
	/*------------------Backward projection------------------------------*/
    	for(i=0;i<prodDims;i++){	temp[i]=0; }
    	backwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, temp, angles, &options, n_angles, n_rays_y, n_rays_z, dims[2], grid_offset, voxel_size);
    	for(i=0;i<prodDims;i++){	Nablafyk[i] = Nablafyk[i] + temp[i]; } /* For each voxel */

	/* Take the projected step from yk to xkp1 */
	t = - 1/bL; /* bL is original setting of L_k */
	dcopy_(&prodDims,yk,&one,xkp1,&one); 
	daxpy_(&prodDims,&t,Nablafyk,&one,xkp1,&one); /* x_k+1 = y_k - t*Nablaf(y_k) */
	
	P(xkp1,ctype,d,c,prodDims);


	/* Backtracking on Lipschitz parameter. */
	for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-yk[i]; }
	
	hxkp1[0] = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D); /* alpha*T_tau(xkp1+1) (Nablafyk is returned as gradient) */
    
	/*-----------------Forward projection--------------------------------*/
	for(i=0;i<n_rays;i++){	tv2[i]=0; }
	forwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, xkp1 , angles, &options, n_angles, n_rays_y, n_rays_z); /*tv2 = A*x_k+1 */
	for(i=0;i<n_rays;i++){ tv2[i] = tv2[i]-b[i]; } /* tv2 = (A*x_k+1 - b) */
	
	gxkp1[0] = 0.5*pow( dnrm2_(&n_rays,tv2,&one),2); /* 0.5*||A*x_k+1 - b||^2 */
	
	numFunc[0]+=1;
	fxkp1[0] = hxkp1[0] + gxkp1[0];		/* f(x_k+1) = h(x_k+1) + g(x_k+1) = alpha*T_tau(x_k+1) + 0.5*||A*x_k+1 - b||^2 */
	
	while( fxkp1[0]/(1+1e-14) > fyk + ddot_(&prodDims,Nablafyk,&one,tv,&one) + (bL/2)*pow(dnrm2_(&prodDims,tv,&one),2) ){
		numBack[0]+=1;
		bL = s_L*bL;
				
		/* Take the projected step from yk to xkp1 */
		t = - 1/bL; /* t = - L^-1 */
		dcopy_(&prodDims,yk,&one,xkp1,&one); /* x_k+1 = y_k */
		daxpy_(&prodDims,&t,Nablafyk,&one,xkp1,&one); /* x_k+1 = x_k+1 - t*Nablafyk */
		P(xkp1,ctype,d,c,prodDims); /* project to within bounds (x_k+1) */
		
		/* Backtracking on Lipschitz parameter. */
		for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-yk[i]; }
			
		hxkp1[0] = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D); /* alpha*T_tau(x_k+1) */
		
		/*-----------------Forward projection--------------------------------*/
		for(i=0;i<n_rays;i++){	tv2[i]=0; }
        	forwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, xkp1, angles, &options, n_angles, n_rays_y, n_rays_z); /*tv2 = A*x_k+1 */
        	for(i=0;i<n_rays;i++){ tv2[i] = tv2[i]-b[i]; } /* tv2 = (A*x_k+1 - b) */
			
		gxkp1[0] = 0.5*pow( dnrm2_(&n_rays,tv2,&one),2); /* 0.5*||A*x_k+1 - b||^2 */
		
		numFunc[0]+=1;
		fxkp1[0] = hxkp1[0] + gxkp1[0]; /* f(x_k+1) = h(x_k+1) + g(x_k+1) = alpha*T_tau(x_k+1) + 0.5*||A*x_k+1 - b||^2 */
	}


	/* Calculate initial gradient map */
	if(ctype==1){ /* No constraints - just real numbers */
		nGt = dnrm2_(&prodDims,Nablafxkp1,&one);
	}
	else{ /* Positivity constraint */
		t = - 1/bL;
		dcopy_(&prodDims,xkp1,&one,tv,&one);
		daxpy_(&prodDims,&t,Nablafxkp1,&one,tv,&one);
		P(tv,ctype,d,c,prodDims);
		
		for(i=0;i<prodDims;i++){	tv[i] = xkp1[i]-tv[i]; }
		
		nGt = bL*dnrm2_(&prodDims,tv,&one);
	}

	/* save the initial parameter */
	Lm1 = bL;
	nGtm1 = nGt;


	dcopy_(&prodDims,xkp1,&one,xk,&one);

	dcopy_(&prodDims,xkp1,&one,yk,&one);
    
	
  /* LOOP */
	stop = 0; /*Flag for when to break the for-loop*/

	/*Calculate fxk */
	fxk = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
		
	/*-----------------Forward projection--------------------------------*/
	for(i=0;i<n_rays;i++){	tv2[i]=0; }
	forwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, xkp1, angles, &options, n_angles, n_rays_y, n_rays_z);
	for(i=0;i<n_rays;i++){ tv2[i] = tv2[i]-b[i]; }

	numFunc[0]+=1;
	fxk = fxk + 0.5*pow( dnrm2_(&n_rays,tv2,&one),2);


  	for(k=0; k<=k_max; k++){
		kk+=1;
		Lklist[kk]=bL;
		muklist[kk]=bmu;
		/* Calculate the gradient in yk */
		numGrad[0]+=1;
		fyk = alpha*DTD(yk,Nablafyk,uijl,tau,D);

		for(i=0;i<prodDims;i++){	Nablafyk[i]=alpha*Nablafyk[i]; }

		/*-----------------Forward projection--------------------------------*/
		for(i=0;i<n_rays;i++){ tv2[i]=0; }
        	forwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, yk, angles, &options, n_angles, n_rays_y, n_rays_z);
		for(i=0;i<n_rays;i++){	tv2[i] = tv2[i]-b[i]; }
		
		numFunc[0]+=1;
		fyk = fyk + 0.5*pow( dnrm2_(&n_rays,tv2,&one),2);

		for(i=0;i<prodDims;i++){	tv[i]=0; }
	
		/*------------------Backward projection------------------------------*/
		for(i=0;i<prodDims;i++){	temp[i]=0; }
	        backwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, temp, angles, &options, n_angles, n_rays_y, n_rays_z, dims[2], grid_offset, voxel_size);
        	for(i=0;i<prodDims;i++){	Nablafyk[i] = Nablafyk[i] + temp[i]; } /* For each voxel */
	
		/* Update estimate of the strong convexity parameter as minimum of current value and the computed value between xk and yk */
		if(k != 0){
			for(i=0;i<prodDims;i++){ tv[i] = xk[i]-yk[i]; }

			bmu = MAX( MIN(2*(fxk*(1+1e-14)-(fyk+ddot_(&prodDims,Nablafyk,&one,tv,&one)))/pow(dnrm2_(&prodDims,tv,&one),2),bmu),0);
		}

		/* Take the projected step from yk to xkp1 */
		t = - 1/bL;
		dcopy_(&prodDims,yk,&one,xkp1,&one);
		daxpy_(&prodDims,&t,Nablafyk,&one,xkp1,&one);

		P(xkp1,ctype,d,c,prodDims);
		
		/* Backtracking on Lipschitz parameter. */
		for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-yk[i]; }
		
		hxkp1[0] = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
		
		/*-----------------Forward projection--------------------------------*/
		for(i=0;i<n_rays;i++){	tv2[i]=0; }
		forwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, xkp1, angles, &options, n_angles, n_rays_y, n_rays_z);
	        for(i=0;i<n_rays;i++){ tv2[i] = tv2[i]-b[i]; }

		gxkp1[0] = 0.5*pow( dnrm2_(&n_rays,tv2,&one),2);
		
		numFunc[0]+=1;
		fxkp1[0] = hxkp1[0] + gxkp1[0];		

		while( fxkp1[0]/(1+1e-14) > fyk + ddot_(&prodDims,Nablafyk,&one,tv,&one) + (bL/2)*pow(dnrm2_(&prodDims,tv,&one),2) ){
			numBack[0]+=1;
			bL = s_L*bL;
				
			/* Take the projected step from yk to xkp1 */
			t = - 1/bL;
			dcopy_(&prodDims,yk,&one,xkp1,&one);
			daxpy_(&prodDims,&t,Nablafyk,&one,xkp1,&one);
			P(xkp1,ctype,d,c,prodDims);
			
			/* Backtracking on Lipschitz parameter. */
			for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-yk[i]; }
			
			hxkp1[0] = alpha*DTD(xkp1,Nablafxkp1,uijl,tau,D);
			
			/*-----------------Forward projection--------------------------------*/
			for(i=0;i<n_rays;i++){	tv2[i]=0; }
        		forwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, xkp1, angles, &options, n_angles, n_rays_y, n_rays_z);
			for(i=0;i<n_rays;i++){ tv2[i] = tv2[i]-b[i]; }
			
			gxkp1[0] = 0.5*pow( dnrm2_(&n_rays,tv2,&one),2);

			numFunc[0]+=1;
			fxkp1[0] = hxkp1[0] + gxkp1[0];
		}
		
		fxkp1l[kk] = fxkp1[0];
		
		if(ghxl){
			hxkp1l[kk] = hxkp1[0];
			gxkp1l[kk] = gxkp1[0];
		}
			
		/* store the iterate if requested */
		if(xl){
			for(i=0;i<prodDims;i++){
				xlist[kk*prodDims+i] = xkp1[i];
			}
		}


		if(verbose)
			printf("k=%6d  f(x^k+1)=%e  ||G_L(x^k+1)||=%e  L_k=%.2e  mu_k=%.2e\n",kk,fxkp1[0],nGt,bL,bmu);DRAW;

		/* calculate the gradient in xkp1 */
		numGrad[0]+=1;
		for(i=0;i<prodDims;i++){	Nablafxkp1[i]=alpha*Nablafxkp1[i]; }
        
		/*------------------Backward projection------------------------------*/
		for(i=0;i<prodDims;i++){	temp[i]=0; }
	        backwardProjection(source_x, source_y, source_z, det_x, det_y, det_z, tv2, temp, angles, &options, n_angles, n_rays_y, n_rays_z, dims[2], grid_offset, voxel_size);
        	for(i=0;i<prodDims;i++){	Nablafxkp1[i] = Nablafxkp1[i] + temp[i]; } /* For each voxel */

		/* Check stopping criteria xkp1*/
		if(ctype==1){
			nGt = dnrm2_(&prodDims,Nablafxkp1,&one);
			if(nGt <= epsb_rel*prodDims){
				stop = 1;
   			/*overwrite xkp1 to return with*/
			  t = - 1/bL;
				daxpy_(&prodDims,&t,Nablafxkp1,&one,xkp1,&one); 
			}
		}
		else{
			t = - 1/bL;
			dcopy_(&prodDims,xkp1,&one,tv,&one);
			daxpy_(&prodDims,&t,Nablafxkp1,&one,tv,&one);
			P(tv,ctype,d,c,prodDims);

			for(i=0;i<prodDims;i++){	tv[i] = xkp1[i]-tv[i]; }
			
			nGt = bL*dnrm2_(&prodDims,tv,&one);
			if(nGt <= epsb_rel*prodDims){
				stop = 1;
   			/*overwrite xkp1 to return with*/
				daxpy_(&prodDims,&t,Nablafxkp1,&one,xkp1,&one);
				P(xkp1,ctype,d,c,prodDims);
			}			
		}

		if(stop || STOPMARK || kk==k_max){
			goto cleanup;
		}

		/* Check stopping criteria yk*/
		if(ctype==1){
			if(dnrm2_(&prodDims,Nablafyk,&one) <= epsb_rel*prodDims)
				stop = 1;
		}
		else{
			t = - 1/bL;
			for(i=0;i<prodDims;i++){	tv[i] = yk[i]-xkp1[i]; }
			if(bL*dnrm2_(&prodDims,tv,&one) <= epsb_rel*prodDims)
				stop = 1;
		}

		if(stop || STOPMARK || kk==k_max){
			goto cleanup;
		}

		/* Compute values for accelerated step size*/
		q = bmu/bL;
		
		if(k!=0)
			cumprod = cumprod*(1-sqrt(q));
		/*tjeck if the convergence rate is fast enough*/
		if(bmu >0){
			if(nGt*nGt> cumprod*(4*bL/bmu-bL/Lm1+4*gamma0*bL/pow(bmu,2))*nGtm1*nGtm1){
			  /*printf("not fast enough %d\n",kk);*/
				/*printf("%f %f %f %f\n",nGt,nGtm1*nGtm1,cumprod,(4*bL/bmu-bL/Lm1+4*gamma0*bL/pow(bmu,2)));*/
				/*list of restart positions */
				p_rp = AddItem(p_rp,kk);
				bmu=bmu*s_mu;
				numRest[0]+=1;
				goto restart;
			}
		}
		/*printf("%f\n",q);DRAW;*/

		thetakp1 = (-(pow(thetak,2)-q)+sqrt(pow( pow(thetak,2)-q,2)+4*pow(thetak,2)))/2.0;
		betak = (thetak*(1-thetak))/(pow(thetak,2)+thetakp1);

		if(k==0){
			gamma0 = thetakp1*(thetakp1*bL-bmu)/(1-thetakp1);
			/*printf("gamma0 %f\n",gamma0);*/
		}

		/* accelerated term*/
		/* yk = xkp1 + betak*(xkp1-xk) */
		for(i=0;i<prodDims;i++){ tv[i] = xkp1[i]-xk[i]; }

		dcopy_(&prodDims,xkp1,&one,yk,&one);
		daxpy_(&prodDims,&betak,tv,&one,yk,&one);

		/* Update the values for the next iteration*/
		thetak = thetakp1;
		dcopy_(&prodDims,xkp1,&one,xk,&one);
        
		fxk = fxkp1[0];

	}

	cleanup:
  	free(Nablafyk);
  	free(Nablafxkp1);

	free(yk);
	free(xk);

	free(tv);
	free(tv2);
	free(uijl);
    	free(temp);

	kend[0]=(double)kk;

}
