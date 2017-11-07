#include "defs.h"


bool fbp_gpu_20_fs(allData *aD){
	unsigned int onx, ony, wo, ho, obnx, obny, nx, ny, wh, wh2, wu, vH;
	float sx, sy, ps;
	double angleStart, angleStep;
	Ipp32f *vux, *vuy, *v;
	Ipp64f *v1, *v2, *v3, *v4;
	xFBP_gpu *fg;
	gpu_info *gi;
	xData *xd;
	xRoi *xr;

	fg = aD->fg;
	gi = aD->gi;
	xd = aD->data;
	xr = aD->hms->fbp->backprojection->roi;

	nx = aD->ta_nx;
	ny = aD->ta_ny;
	wh = gi->ho * gi->wo;
	wh2 = wh/2;
	wu = ny/wh2;
	
	onx = gi->mxo;
	ony = gi->myo;
	wo = gi->wo;
	ho = gi->ho;
	obnx = onx/wo;
	obny = ony/ho;

	ps = gi->outputPixelSize;

	angleStart = pi*double(xr->angle)/180.0;
	angleStep = double(gi->rotAngleStep); 
	
	v1 = ippsMalloc_64f(ny);
	v2 = ippsMalloc_64f(ny);
	v3 = ippsMalloc_64f(ny);
	v4 = ippsMalloc_64f(2*ny);
		
	xd->vecCS = ippsMalloc_32f(2*ny);
	ippsVectorSlope_64f(v1,ny,angleStart,angleStep);
	ippsCos_64f_A53(v1,v2,ny);
	ippsSin_64f_A53(v1,v3,ny);
	
	ippsRealToCplx_64f(v2, v3, (Ipp64fc *)v4, ny);
	ippsConvert_64f32f(v4,xd->vecCS,2*ny);

	ippsFree(v1); v1 = NULL;
	ippsFree(v2); v2 = NULL;
	ippsFree(v3); v3 = NULL;
	ippsFree(v4); v4 = NULL;



	sx = 0.5f*ps * float(onx - 1);
	sy = 0.5f*ps * float(ony - 1);

	aD->roiRadius = sqrt(sx*sx+sy*sy);

	aD->to_nx = onx;
	aD->to_ny = ony;

	vux = ippsMalloc_32f(onx);
	vuy = ippsMalloc_32f(onx);

	if(vux == NULL || vuy == NULL){
		printError(aD,"can not allocate memory (fbp, ux, uy)");
		return false;
	}

	xd->vecXY_block = ippsMalloc_32f(2*wo*ho);
	xd->vecXY = ippsMalloc_32f(2*onx*ony);
	xd->vecRA = ippsMalloc_32f(2*onx*ony);
	xd->vecbXY = ippsMalloc_32f(2*obnx*obny);
		
	if(xd->vecXY == NULL || xd->vecbXY == NULL || xd->vecRA == NULL){
		printError(aD, "can not allocate memory (FBP, vecXY, vecbXY, vecRA)");
		return false;
	}

	ippsVectorSlope_32f(vux,wo,-0.5f*ps*(wo-1),ps);
	for(unsigned int i=0; i<ho; i++){
		ippsSet_32f(-0.5f*ps*(ho-1)+float(i)*ps,vuy,wo);
		ippsRealToCplx_32f(vux,vuy,(Ipp32fc*)(xd->vecXY_block)+i*wo,wo);
	}

	ippsVectorSlope_32f(vux,onx,-0.5f*ps*(onx-1),ps);
	for(unsigned int i=0; i<ony; i++){
		ippsSet_32f(-0.5f*ps*(ony-1)+float(i)*ps,vuy,onx);
		ippsRealToCplx_32f(vux,vuy,(Ipp32fc*)(xd->vecXY)+i*onx,onx);
	}

	ippsCartToPolar_32fc((Ipp32fc *)xd->vecXY, xd->vecRA, xd->vecRA+onx*ony, onx*ony);
	v =  xd->vecRA+onx*ony;
	for(unsigned int i=0; i<onx*ony; i++){
		if(v[i]<0.0f) v[i] += (2.0f*pif);
	}
	ippsMulC_32f_I(float(aD->ta_ny_13-1)/pif,v,onx*ony);
	ippsThreshold_LTValGTVal_32f_I(v, onx*ony, 0.0f, 0.0f, float(2*aD->ta_ny_13-2), float(2*aD->ta_ny_13-2));
			
	ippsFree(vux); vux = NULL;
	ippsFree(vuy); vuy = NULL;

	xd->vecAT = ippsMalloc_32f(onx*ony);
	ippsPhase_32fc((Ipp32fc*)xd->vecXY, xd->vecAT, onx*ony);
	ippsAddC_32f_I(pif,xd->vecAT,onx*ony);
	ippsMulC_32f_I(float(aD->ta_ny_13)/pif,xd->vecAT,onx*ony);


	for(unsigned int i=0; i<obny; i++){
		for(unsigned int j=0; j<obnx; j++){
			xd->vecbXY[2*(i*obnx+j)] = 0.5f * (xd->vecXY[2*(i*ho*onx+j*wo)] + xd->vecXY[2*(((i+1)*ho-1)*onx+(j+1)*wo-1)]);
			xd->vecbXY[2*(i*obnx+j)+1] = 0.5f * (xd->vecXY[2*(i*ho*onx+j*wo)+1] + xd->vecXY[2*(((i+1)*ho-1)*onx+(j+1)*wo-1)+1]);
		}
	}
		
	xd->vto = ippsMalloc_32f(onx*ony);
	if(xd->vto == NULL){
		printError(aD, "can not allocate memory (output matrix)");
		return false;
	}

	if(!form_filter(aD))return false;

	vH = gi->vertH;
	xd->veccF = ippsMalloc_32fc(nx*vH);
	
	Ipp32f *vz = ippsMalloc_32f(nx);
	ippsZero_32f(vz, nx);
		
	ippsRealToCplx_32f(xd->vfilter, vz, xd->veccF, nx);
	for(unsigned int i=1; i<vH; i++){
		ippsCopy_32fc(xd->veccF, xd->veccF + i*nx, nx);
	}
	ippsFree(vz); vz = NULL;

	
	return true;
}


bool fbp_gpu_fs(allData *aD){
	xData *xd;
	xFBP_gpu *fg;
	gpu_info *gi;
	xRoi *xr;

	unsigned int onx, ony, wo, ho, ny, nx;
	double angleStart, angleStep;
	unsigned int ont;
	unsigned int ndiam, ndiameter, nleft;
	unsigned int vH;
	float sx, sy, ps;
	double dps, xmin, ymin, yval;
	double rat;
	Ipp64f *v1, *v2, *v3, *v4;

	xd = aD->data;
	fg = aD->fg;
	gi = aD->gi;

	fg->nx = aD->ta_nx;
	fg->ny = aD->ta_ny;
	fg->nxo = gi->mxo;
	fg->nyo = gi->myo;
	fg->blockWidth = gi->wo;
	fg->blockHeight = gi->ho;
	fg->vH = gi->vertH;

	nx = aD->ta_nx;
	ny = aD->ta_ny;
	
	xr = aD->hms->fbp->backprojection->roi;

	ps = gi->outputPixelSize;
	
	angleStart = pi*double(xr->angle)/180.0;
	angleStep = double(gi->rotAngleStep); 
	
	v1 = ippsMalloc_64f(ny);
	v2 = ippsMalloc_64f(ny);
	v3 = ippsMalloc_64f(ny);
	v4 = ippsMalloc_64f(2*ny);
		
	xd->vecCS = ippsMalloc_32f(2*ny);
	ippsVectorSlope_64f(v1,ny,angleStart,angleStep);
	ippsCos_64f_A53(v1,v2,ny);
	ippsSin_64f_A53(v1,v3,ny);
	
	ippsRealToCplx_64f(v2, v3, (Ipp64fc *)v4, ny);
	ippsConvert_64f32f(v4,xd->vecCS,2*ny);

	ippsFree(v1); v1 = NULL;
	ippsFree(v2); v2 = NULL;
	ippsFree(v3); v3 = NULL;
	ippsFree(v4); v4 = NULL;


	onx = gi->mxo;
	ony = gi->myo;
	ont = onx*ony;

	wo = gi->wo;
	ho = gi->ho;

	sx = 0.5f*ps * float(onx - 1);
	sy = 0.5f*ps * float(ony - 1);

	aD->roiRadius = sqrt(sx*sx+sy*sy);

	ndiam = __min(nx,(unsigned int)(ceil(2.0f*aD->roiRadius))+2);
	int_ratio(ndiam,wo,&ndiameter);
	fg->chunkWidth =  ndiameter;
	nleft = (unsigned int)(floor(aD->new_xc - 0.5f*float(ndiameter) - 1.0f));
	if (nleft + fg->chunkWidth > nx){
		nleft = nx - fg->chunkWidth;
	}
	fg->chunkLeft = nleft;
	fg->xc = aD->new_xc - float(nleft);

	fg->numChunks = fg->ny/fg->vH;

	aD->to_nx = onx;
	aD->to_ny = ony;

	xd->vecX = ippsMalloc_32f(ont);
	xd->vecY = ippsMalloc_32f(ont);
	
	dps = double(ps);
	xmin = -0.5 * double(onx-1) * dps;
	ymin = -0.5 * double(ony-1) * dps;

	v1 = ippsMalloc_64f(onx);
	v2 = ippsMalloc_64f(onx);
	ippsVectorSlope_64f(v1,onx,xmin,dps);
	
	rat = double(aD->ta_ny_13-1)/pi;

	for(unsigned int i = 0; i<ony; i++){
		yval = ymin + double(i)*dps;
		ippsSet_64f(yval,v2,onx);
		ippsConvert_64f32f(v1,xd->vecX+i*onx,onx);
		ippsConvert_64f32f(v2,xd->vecY+i*onx,onx);
	}

	ippsFree(v1); v1 = NULL;
	ippsFree(v2); v2 = NULL;
			
	xd->vto = ippsMalloc_32f(onx*ony);
	if(xd->vto == NULL){
		printError(aD, "can not allocate memory (output matrix)");
		return false;
	}
	
	vH = gi->vertH;

	if(!form_filter(aD))return false;
	xd->veccF = ippsMalloc_32fc(nx*vH);
	
	Ipp32f *vz = ippsMalloc_32f(nx);
	ippsZero_32f(vz, nx);
		
	ippsRealToCplx_32f(xd->vfilter, vz, xd->veccF, nx);
	for(unsigned int i=1; i<vH; i++){
		ippsCopy_32fc(xd->veccF, xd->veccF + i*nx, nx);
	}
	ippsFree(vz); vz = NULL;

	return true;
}


bool fbp_gpu_cpu2_fs(allData *aD){
	xData *xd;
	xFBP_gpu *fg;
	gpu_info *gi;
	xRoi *xr;

	unsigned int onx, ony, wo, ho, ny, nx;
	double angleStart, angleStep;
	unsigned int ont;
	unsigned int ndiam, ndiameter, nleft;
	unsigned int vH;
	float sx, sy, ps;
	double dps, xmin, ymin, yval;
	double rat;
	Ipp64f *v1, *v2, *v3, *v4;

	xd = aD->data;
	fg = aD->fg;
	gi = aD->gi;

	fg->nx = aD->ta_nx;
	fg->ny = aD->ta_ny;
	fg->nxo = gi->mxo;
	fg->nyo = gi->myo;
	fg->blockWidth = gi->wo;
	fg->blockHeight = gi->ho;
	fg->vH = gi->vertH;

	nx = aD->ta_nx;
	ny = aD->ta_ny;
	
	xr = aD->hms->fbp->backprojection->roi;

	ps = gi->outputPixelSize;

	unsigned int na, nr, pol_a, pol_r;

	wo = gi->wo;
	ho = gi->ho;

	na = 2*aD->ta_ny_13 -1;
	int_ratio(na,wo,&pol_a);

	
	angleStart = pi*double(xr->angle)/180.0;
	angleStep = double(gi->rotAngleStep); 
	
	v1 = ippsMalloc_64f(ny);
	v2 = ippsMalloc_64f(ny);
			
	xd->vecCos = ippsMalloc_32f(aD->ta_ny_13);
	ippsVectorSlope_64f(v1,ny,angleStart,angleStep);
	ippsCos_64f_A53(v1,v2,aD->ta_ny_13);
	
	ippsConvert_64f32f(v2,xd->vecCos,aD->ta_ny_13);

	ippsFree(v1); v1 = NULL;
	ippsFree(v2); v2 = NULL;

	onx = gi->mxo;
	ony = gi->myo;
	
	sx = 0.5f*ps * float(onx - 1);
	sy = 0.5f*ps * float(ony - 1);

	aD->roiRadius = sqrt(sx*sx+sy*sy);

	nr = (unsigned int)(ceil(aD->roiRadius/ps)+1.0f);

	int_ratio(nr,ho,&pol_r);

	xd->vecPol = ippsMalloc_32f(pol_r*pol_a);

	gi->pol_r = pol_r;
	gi->pol_a = pol_a;

	onx = gi->mxo;
	ony = gi->myo;
			
	xd->vto = ippsMalloc_32f(onx*ony);
	if(xd->vto == NULL){
		printError(aD, "can not allocate memory (output matrix)");
		return false;
	}
	
	vH = gi->vertH;

	if(!form_filter(aD))return false;
	xd->veccF = ippsMalloc_32fc(nx*vH);
	
	Ipp32f *vz = ippsMalloc_32f(nx);
	ippsZero_32f(vz, nx);
		
	ippsRealToCplx_32f(xd->vfilter, vz, xd->veccF, nx);
	for(unsigned int i=1; i<vH; i++){
		ippsCopy_32fc(xd->veccF, xd->veccF + i*nx, nx);
	}
	ippsFree(vz); vz = NULL;

	return true;
}

bool fbp_gpu(allData *aD){
	printTagStart(aD,"FilteredBackprojectionGPU");

#ifdef _ALGORITHM2
	if(aD->isFirstSlice){
		if(!fbp_gpu_20_fs(aD))return false;
	}
	if(!fbp_cuda_20(aD))return false;
#elif _ALGORITHM_CPU_2
	if(aD->isFirstSlice){
		if(!fbp_gpu_cpu2_fs(aD))return false;
	}
	if(!fbp_cuda_cpu2(aD))return false;

#else
	if(aD->isFirstSlice){
		if(!fbp_gpu_fs(aD))return false;
	}
	if(!fbp_cuda(aD))return false;
#endif

	printTagEnd(aD);
	return true;
}







