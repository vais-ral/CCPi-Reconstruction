
#include <mex.h>
#include <list>
#include "mex_types.hpp"
#include "instruments.hpp"
#include "algorithms.hpp"
#include "blas.hpp"

#if (defined(_WIN32) || defined(__WIN32__) )
#define DRAW mexEvalString("drawnow;");
#else
#define DRAW ;
#endif

struct Dtype {
        int dim;
        int m,n,l;
        long prodDims;
};

static void dcopyf(const long n, const float x[], const int incx, real y[],
		   const int incy);

// dummy cone-beam device
namespace CCPi {

  class dummy_cone : public cone_beam {
  public:
    bool setup_experimental_geometry(const std::string path,
				     const std::string file,
				     const bool phantom = false);
    bool read_scans(const std::string path, const bool phantom = false);
    bool finish_voxel_geometry(real voxel_origin[3], real voxel_size[3],
			       const voxel_data &voxels) const;
    void apply_beam_hardening();
  };
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  register double *b,*d,*c,*dims,alpha,tau,bL,bmu,epsb_rel,*hxkp1l,*gxkp1l,*xlist, *voxel_size;
  register double *source_x, *source_y, *source_z, *det_x, *det_y, *det_z, *angles, *grid_offset; 
  register double *fxkp1,*hxkp1,*gxkp1,*fxkp1l,*k,*numGrad,*numBack,*numFunc,*numRest,*Lklist,*muklist,*rklist;
  register float *x, *xkp1;
  mxArray *M,*S,*Mdims;
  int i,j,k_max,dim,ctype,ghxl,xl,verbose,temp, n_rays_y, n_rays_z, n_angles;
  long prodDims;
  Dtype D;

  if(nrhs != 24){
    printf("Should contain 24 input parameters but has %i\n",nrhs); DRAW;}
  else{							

    std::list<int> rp;

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

    dim = std::max( mxGetM(M), mxGetN(M) );
                
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
    dim = std::max( mxGetM(Mdims), mxGetN(Mdims) );
    prodDims=1;
    for(i=0;i<dim;i++){
      prodDims = prodDims*(int)dims[i];
    }
    D.dim = dim;
    D.m=(int)dims[0];D.n=(int)dims[1];
    if(dim==3){
      D.l=(int)dims[2];
    }else
      D.l = 1;
    D.prodDims =prodDims;
        
    /*Allocate memory and assign output pointer*/
    plhs[0] = mxCreateNumericMatrix(prodDims, 1, mxSINGLE_CLASS, mxREAL); 
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
    xkp1 = (float*)mxGetPr(plhs[0]);
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

    dcopy(prodDims,x,1,xkp1,1); 

    CCPi::cone_beam *dev = new CCPi::dummy_cone;
    dev->set_params(*source_x, *source_y, *source_z, *det_x, det_y, det_z,
		    angles, n_rays_y, n_rays_z, n_angles);

    int ki = 0;
    CCPi::tvreg_core(xkp1,fxkp1,hxkp1,gxkp1,fxkp1l,&ki,voxel_size, b,alpha,
		     tau,bL,bmu,epsb_rel,k_max,D.dim, D.m, D.n, D.l, D.prodDims,
		     ctype,d,c,(bool)ghxl,(bool)xl,hxkp1l,gxkp1l,xlist,
		     (bool)verbose,numGrad,numBack,numFunc,numRest,Lklist,
		     muklist,rp,grid_offset, dev);

    delete dev;
    *k = (double)ki;
    /*write the dynamical allocated restart list to a vector with the correct dimensions*/			
    long rql = (long)rp.size();

    plhs[15] = mxCreateDoubleMatrix( rql-1, 1, mxREAL);
    rklist = mxGetPr(plhs[15]);

    rql = 0;
    for (std::list<int>::const_iterator p = rp.begin(); p != rp.end(); ++p) {
      rklist[rql] = *p;
      rql++;
    }
  }

}

void dcopyf(const long n, const float x[], const int incx, real y[],
	   const int incy)
{
  for (long i = 0; i < n; i++)
    y[i] = real(x[i]);
}

bool CCPi::dummy_cone::setup_experimental_geometry(const std::string path,
						   const std::string file,
						   const bool phantom)
{
  return false;
}

bool CCPi::dummy_cone::read_scans(const std::string path, const bool phantom)
{
  return false;
}

bool CCPi::dummy_cone::finish_voxel_geometry(real voxel_origin[3],
					     real voxel_size[3],
					     const voxel_data &voxels) const
{
  return false;
}

void CCPi::dummy_cone::apply_beam_hardening()
{
}
