#include <mex.h>
#include "tools.h"
#include "tv_core.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	register double *b,*d,*c,*dims,alpha,tau,bL,bmu,epsb_rel,*hxkp1l,*gxkp1l,*xlist, *voxel_size;
    register double *source_x, *source_y, *source_z, *det_x, *det_y, *det_z, *angles, *grid_offset; 
	register double *xkp1,*fxkp1,*hxkp1,*gxkp1,*fxkp1l,*k,*numGrad,*numBack,*numFunc,*numRest,*Lklist,*muklist,*rklist;
    register float *x;
	mxArray *M,*S,*Mdims;
	int i,j,k_max,dim,prodDims,one=1,ctype,ghxl,xl,verbose,rql,temp, n_rays_y, n_rays_z, n_angles;  
	Dtype D;
	listelement *p_rp;

	p_rp = NULL;
	p_rp = AddItem(p_rp,0);  /*Need to initialize to avoid pointer problems when passing to the core function*/
    
	if(nrhs != 24){
		printf("Should contain 24 input parameters but has %i\n",nrhs); DRAW;}
  else{							

        Mdims = (mxArray*)prhs[0];
		voxel_size = mxGetPr(Mdims);
	
		M = (mxArray*)prhs[1];
		b = mxGetPr(M);

		S = (mxArray*)prhs[2];
		alpha = (double)(mxGetScalar(S));
		
		S = (mxArray*)prhs[3];
		tau = (double)(mxGetScalar(S));
		
		Mdims = (mxArray*)prhs[4];
		dims = mxGetPr(Mdims); 
		
		S = (mxArray*)prhs[5];
		bL = (double)(mxGetScalar(S));
		
		S = (mxArray*)prhs[6];
		bmu = (double)(mxGetScalar(S));
		
		S = (mxArray*)prhs[7];
		epsb_rel = (double)(mxGetScalar(S));

		S = (mxArray*)prhs[8];
		k_max = (int)(mxGetScalar(S));

		M = (mxArray*)prhs[9];
		x = (float*)mxGetData(M);

        dim = MAX( mxGetM(M), mxGetN(M) );
                
		S = (mxArray*)prhs[10];
		ctype = (int)(mxGetScalar(S));

		M = (mxArray*)prhs[11];
		d = mxGetPr(M);

		M = (mxArray*)prhs[12];
		c = mxGetPr(M);

		S = (mxArray*)prhs[13];
		ghxl = (int)(mxGetScalar(S));

		S = (mxArray*)prhs[14];
		xl = (int)(mxGetScalar(S));

		S = (mxArray*)prhs[15];
		verbose = (int)(mxGetScalar(S));
            
     	source_x = mxGetPr(prhs[16]); 
        source_y = mxGetPr(prhs[17]); 
        source_z = mxGetPr(prhs[18]); 
        det_x = mxGetPr(prhs[19]); 
        det_y = mxGetPr(prhs[20]); 
        det_z = mxGetPr(prhs[21]); 
        angles = mxGetPr(prhs[22]);
        grid_offset = mxGetPr(prhs[23]);
        
        n_rays_y = mxGetM(prhs[20]);
        n_rays_z = mxGetM(prhs[21]);
        n_angles = mxGetM(prhs[22]);
        
        /*obtain the dimensions */
		dim = MAX( mxGetM(Mdims), mxGetN(Mdims) );
		prodDims=1;
		for(i=0;i<dim;i++){
			prodDims = prodDims*(int)dims[i];
        }
		D.dim = dim;
		D.m=(int)dims[0];D.n=(int)dims[1];
		if(dim==3){
			D.l=(int)dims[2];
		}
		D.prodDims =prodDims;
        
		/*Allocate memory and assign output pointer*/
		plhs[0] = mxCreateDoubleMatrix(prodDims, 1, mxREAL); 
		plhs[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
		plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
		plhs[3] = mxCreateDoubleMatrix(1, 1, mxREAL);
		plhs[4] = mxCreateDoubleMatrix(k_max+1, 1, mxREAL);
    	plhs[5] = mxCreateDoubleMatrix(1, 1, mxREAL);
			
		if(ghxl){
			plhs[6] = mxCreateDoubleMatrix(k_max+1, 1, mxREAL);
			plhs[7] = mxCreateDoubleMatrix(k_max+1, 1, mxREAL);}
		else{
				plhs[6] = mxCreateDoubleMatrix(1, 1, mxREAL);
				plhs[7] = mxCreateDoubleMatrix(1, 1, mxREAL);
			}

			if(xl){
				plhs[8] = mxCreateDoubleMatrix(prodDims*(k_max+1), 1, mxREAL);
			}
			else{
				plhs[8] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			}

			plhs[9] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			plhs[10] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			plhs[11] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			plhs[12] = mxCreateDoubleMatrix( 1, 1, mxREAL);
			plhs[13] = mxCreateDoubleMatrix( k_max+1, 1, mxREAL);
			plhs[14] = mxCreateDoubleMatrix( k_max+1, 1, mxREAL);

			/* Get a pointer to the data space in our newly allocated memory */
			xkp1 = mxGetPr(plhs[0]);
			fxkp1 = mxGetPr(plhs[1]);
			hxkp1 = mxGetPr(plhs[2]);
			gxkp1 = mxGetPr(plhs[3]);
			fxkp1l = mxGetPr(plhs[4]);
			k = mxGetPr(plhs[5]);
			hxkp1l = mxGetPr(plhs[6]);
			gxkp1l = mxGetPr(plhs[7]);
			xlist = mxGetPr(plhs[8]);

			numGrad = mxGetPr(plhs[9]);
			numBack = mxGetPr(plhs[10]);
			numFunc = mxGetPr(plhs[11]);
			numRest = mxGetPr(plhs[12]);
			Lklist = mxGetPr(plhs[13]);
			muklist = mxGetPr(plhs[14]);

            dcopyf_(&prodDims,x,&one,xkp1,&one); 
			
			tvreg_upn_core(xkp1,fxkp1,hxkp1,gxkp1,fxkp1l,k,voxel_size, dims,b,alpha,tau,bL,bmu,epsb_rel,k_max,D,ctype,d,c,
                           ghxl,xl,hxkp1l,gxkp1l,xlist,verbose, numGrad,numBack,numFunc,numRest,Lklist,muklist,p_rp,
                           source_x, source_y, source_z, det_x, det_y, det_z, angles, grid_offset, n_rays_y, n_rays_z, n_angles);

			/*write the dynamical allocated restart list to a vector with the correct dimensions*/			
			rql = QueueLength(p_rp);


			plhs[15] = mxCreateDoubleMatrix( rql-1, 1, mxREAL);
			rklist = mxGetPr(plhs[15]);

			WriteQueueData(p_rp,rklist,rql-1); /*write the list, minus the intial*/
			ClearQueue(p_rp);
	}

}
