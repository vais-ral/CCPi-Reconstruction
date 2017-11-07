/* ----------------------------------------------------------------------
 * PolarToCart.cpp()
 * ----------------------------------------------------------------------
 * Convert between coordinate systems.
 *
 * HMXIF build only. 
 * ----------------------------------------------------------------------
 */

#include "defs.h"

PolarToCart::PolarToCart(){
	vec_x = NULL;
	vec_y = NULL;
	vec_a = NULL;
	vec_r = NULL;
	i_a = NULL;
	i_r = NULL;
	pc_r = NULL;
	pc_a = NULL;
	pc_pos = NULL;
	pc_size = NULL;
	pc_res = NULL;
}

PolarToCart::~PolarToCart(){
	if(vec_x != NULL){ippsFree(vec_x); vec_x = NULL;}
	if(vec_y != NULL){ippsFree(vec_y); vec_y = NULL;}
	if(vec_a != NULL){ippsFree(vec_a); vec_a = NULL;}
	if(vec_r != NULL){ippsFree(vec_r); vec_r = NULL;}
	if(i_a != NULL){ippsFree(i_a); i_a = NULL;}
	if(i_r != NULL){ippsFree(i_r); i_r = NULL;}
	if(pc_r != NULL){ippsFree(pc_r); pc_r = NULL;}
	if(pc_a != NULL){ippsFree(pc_a); pc_a = NULL;}
	if(pc_size != NULL){ippsFree(pc_size); pc_size = NULL;}
	if(pc_pos != NULL){ippsFree(pc_pos); pc_pos = NULL;}
	if(pc_res != NULL){ippsFree(pc_res); pc_res = NULL;}
}

bool PolarToCart::formCudaMatrix(allData *aD){
	xAllFDK *af = aD->hms->af;
	int nx = aD->nx;
	float dd;
	
	dd = (xmax-xmin)/float(OutputWidth-1);
	OutputHeight = int ((ymax-ymin)/dd)  + 1;
		
	int ot = OutputWidth * OutputHeight;
	
	vec_x = ippsMalloc_32f(ot);
	vec_y = ippsMalloc_32f(ot);
	vec_a = ippsMalloc_32f(ot);
	vec_r = ippsMalloc_32f(ot);

	ippsVectorSlope_32f(vec_x,OutputWidth,xmin,dd);

	
	for(int j=1;j<OutputHeight;j++){
		ippsCopy_32f(vec_x, vec_x+j*OutputWidth, OutputWidth);
	}
	for(int i=0;i<OutputHeight;i++){
		ippsSet_32f(ymin + float(i)*dd,vec_y+i*OutputWidth,OutputWidth);
	}

	float xc = float(aD->nx-1)*0.5f + af->axisXShift;

	fprintf(aD->fp,"Image centre is at (%f, centre of the projection) or (%f, absolute value).\n",af->axisXShift,xc);

	ippsSubC_32f_I(xc,vec_x,ot);
	ippsSubC_32f_I(xc,vec_y,ot);
	
	ippsCartToPolar_32f(vec_x, vec_y, vec_r, vec_a, ot);
	ippsFree(vec_x); vec_x=NULL;
	ippsFree(vec_y); vec_y=NULL;


	float rangle = RotAngle;
	
	while(rangle<0.0){
		rangle+=360.0;
	}

	while(rangle>360.0){
		rangle-=360.0;
	}
	rangle*=float(pi/180.0);

	ippsAddC_32f_I(rangle,vec_a,ot);

	float pi2 = (float)(2.0*pi);
	for(int i=0;i<ot;i++){
		if(vec_a[i]<0.0)vec_a[i]+=pi2;
		if(vec_a[i]>=pi2)vec_a[i]-=pi2;
	}

	rmax = __min(xc,float(nx-1-xc));
	
	rmin = 0.0f;
	rmax = __min(rmax,sqrt(__max(xmin*xmin,xmax*xmax) + __max(ymin*ymin,ymax*ymax)));

	if(xmin*xmax>0.0 || ymin*ymax > 0.0){
		ippsMin_32f(vec_r,ot,&rmin);
		rmin-=0.001f;
		if(rmin<=0.0)rmin =0.0;
	}

	ippsThreshold_GTVal_32f_I(vec_r, ot, rmax, -2.0);

	rmax+=0.001f;

	pc_nc = -1;
	bool pq=false;

	float ap = 1.0;

	pc_pos = ippsMalloc_32s(ot);
	pc_size = ippsMalloc_32s(ot);

	for(int i =0;i<ot;i++){
		if(vec_r[i]<-1.0){
			if(!pq)continue;	
			pq = false;
			
		}else{
			if(pq){
				pc_size[pc_nc]++;
			}else{
				pc_nc++;
				pc_size[pc_nc]=1;
				pc_pos[pc_nc]=i;
				pq=true;
			}
		}
	}
	
	pc_nc++;
	printf("pc_nc: %i\n", pc_nc);

	rmin *= ap;
	rmax *= ap;

	printf("rmin: %f\n", rmin);
	printf("rmax: %f\n", rmax);
	
	ippsSum_32s_Sfs(pc_size, pc_nc, &pc_len, 0);

	printf("%i\n", (int)pc_len);

	pc_r = ippsMalloc_32f(pc_len);
	pc_a = ippsMalloc_32f(pc_len);

	int sp = 0;
	for(int i=0;i<pc_nc;i++){
		ippsMulC_32f(vec_r+pc_pos[i],ap,pc_r+sp,pc_size[i]);
		ippsCopy_32f(vec_a+pc_pos[i],pc_a+sp,pc_size[i]);
		sp+=pc_size[i];
	}

	ippsFree(vec_a); vec_a=NULL;
	ippsFree(vec_r); vec_r=NULL;

	Ipp32f amin, amax;

	ippsMinMax_32f(pc_r,pc_len,&amin, &amax);
	printf("pc_r, min: %f, max: %f\n",amin,amax);

	ippsMinMax_32f(pc_a,pc_len,&amin, &amax);

	printf("%f\n",amin);
	printf("%f\n",amax);
	printf("%f\n",0.5*pi);

	int aimax;
	
	aimin = (int) floor(amin*float(num_images)/pi2);
	aimax = (int) ceil(amax*float(num_images)/pi2);
	ailen = aimax-aimin+1;

	int nr = NumberOfCircles;
	int na = ailen;
	
	l_a = w_a *b_a;
	l_r = w_r *b_r;

	nrn = cudaRound(nr,l_r);
	nan = cudaRound(na,l_a);

	printf("l_a: %i\n",l_a);
	printf("l_r: %i\n",l_r);
	printf("na: %i, nan: %i\n",na,nan);
	printf("nr: %i, nrn: %i\n",nr,nrn);

	int man, mrn;
	man = nan/l_a;
	mrn = nrn/l_r;

	Ipp32s *i_w1 = ippsMalloc_32s(pc_len);
	Ipp32s *i_w2 = ippsMalloc_32s(pc_len);
	
	ippsSubC_32f_I(rmin,pc_r,pc_len);

	float dr = (rmax-rmin)/float(NumberOfCircles-1);
	float da = float(2.0*pi)/float(na-1);
	
	ippsDivC_32f_I(dr*float(l_r),pc_r,pc_len);

	ippsConvert_32f32s_Sfs(pc_r, i_w1, pc_len, ippRndZero, 0);
	ippsThreshold_LTValGTVal_32s_I(i_w1, pc_len, 0, 0, mrn-1, mrn-1);

	ippsSubC_32f_I(float(aimin)*da,pc_a,pc_len);
	ippsDivC_32f_I(da*float(l_a),pc_a,pc_len);

	ippsConvert_32f32s_Sfs(pc_a, i_w2, pc_len, ippRndZero, 0);
	ippsThreshold_LTValGTVal_32s_I(i_w2, pc_len, 0, 0, man-1, man-1);

	ippsMulC_32f_I(float(l_r),pc_r,pc_len);
	ippsMulC_32f_I(float(l_a),pc_a,pc_len);
	ippsThreshold_LTValGTVal_32f_I(pc_a, pc_len, 0, 0, float(nan-3), float(nan-3));
	ippsThreshold_LTValGTVal_32f_I(pc_r, pc_len, 0, 0, float(nrn-3), float(nrn-3));

	int mar = man*mrn;
	Ipp32s *i_t = ippsMalloc_32s(mar);

	ippsSet_32s(0,i_t,mar);

	for(int i=0;i<pc_len;i++){
		i_t[i_w1[i]*man+i_w2[i]] = 1;
	}
	ippsFree(i_w1); i_w1 = NULL;
	ippsFree(i_w2); i_w2 = NULL;

	
	ippsSum_32s_Sfs(i_t, mar, &ntb, 0);

	i_a = ippsMalloc_32s(ntb);
	i_r = ippsMalloc_32s(ntb);

	int ij=0;
	for(int i=0;i<mar;i++){
		if(i_t[i] == 0) continue;
		i_a[ij] = i%man;
		i_r[ij] = i/man;
		ij++;
	}
	ippsFree(i_t); i_t = NULL;

	pc_res = ippsMalloc_32f(pc_len);
	ippsZero_32f(pc_res,pc_len);


	return true;
}

