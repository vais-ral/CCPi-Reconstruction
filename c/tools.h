#ifndef __TOOLS__
#define __TOOLS__

#include <mex.h>

#if (defined(_WIN32) || defined(__WIN32__) )
#define DRAW mexEvalString("drawnow;");
#else
#define DRAW ;
#endif

#define MAX(A,B) ((A)>(B)?(A):(B))
#define MIN(A,B) ((A)<(B)?(A):(B))


typedef struct{
	int dim;
	int m,n,l;
	int prodDims;
} Dtype;


typedef struct{
    int     dataitem;
    struct listelement *link;
} listelement;


double min_dbl(double a, double b);

double max_dbl(double a, double b); 

double min3_dbl(double a, double b, double c);

double max3_dbl(double a, double b, double c);

double ceil_j( double arg );

double floor_j( double arg );

int equal_to_precision(double x, double y, double precision);

double alpha_fn(int n, double p1, double p2, double b, double d);

double p(double alpha, double p1, double p2);

double phi(double alpha, double p1, double p2, double b, double d);

double minf(double *x,int N);

double maxf(double *x,int n,int N);

double ddot_(int *n, double *x, int *incx, double *y, int *incy);

double dnrm2_(int *n, double *x, int *incx);

void daxpy_(int *n, double *alpha, double *x, int *incx, double *y, int *incy);

void dcopy_(int *p_n,double *p_x,int *incx,double *p_y,int *incy);

void dcopyf_(int *p_n,float *p_x,int *incx,double *p_y,int *incy);

void P(double *y,int c,double *l,double *u,int N23);

double DTD(double *x,double *Nablafx, double *uijl, double tau, Dtype D);

listelement *AddItem(listelement *listpointer, int data);

void ClearQueue(listelement *listpointer);

void PrintQueue(listelement * listpointer);

int QueueLength(listelement *listpointer);

void WriteQueueData(listelement *listpointer, double *rp, int l);

listelement *RemoveItem(listelement *listpointer);

#endif
