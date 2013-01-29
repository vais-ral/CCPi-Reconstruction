#ifndef __TV_CORE__
#define __TV_CORE__

void tvreg_upn_core(double *xkp1,double *fxkp1,double *hxkp1,double *gxkp1,double *fxkp1l,double *kend, double *voxel_size, double *dims, double *b,double alpha,double tau,
                    double bL,double bmu,double epsb_rel,int k_max,Dtype D,int ctype,double *d,double *c, int ghxl, int xl,double *hxkp1l,double *gxkp1l,double *xlist,
                    int verbose,double *numGrad,double* numBack,double *numFunc,double *numRest,double *Lklist,double *muklist,listelement *p_rp,
                    double *source_x, double *source_y, double *source_z, double *det_x, double *det_y, double *det_z, double *angles, double *grid_offset,
                    int n_rays_y, int n_rays_z, int n_angles);


#endif
