/* ----------------------------------------------------------------------
 * MatPar.cpp()
 * ----------------------------------------------------------------------
 * 
 * ----------------------------------------------------------------------
 */

#include "defs.h"

MatPar::MatPar(){
	m_x = NULL;
	m_y = NULL;
	ivol = NULL;
	rec = NULL;
	vsect = NULL;
	fsect = NULL;
}

MatPar::~MatPar(){
	if(m_x != NULL){ippsFree(m_x); m_x = NULL;}
	if(m_y != NULL){ippsFree(m_y); m_y = NULL;}
	if(ivol != NULL){ippsFree(ivol); ivol = NULL;}
	if(rec != NULL){ippsFree(rec); rec = NULL;}
	if(vsect != NULL){ippsFree(vsect); vsect = NULL;}
	if(fsect != NULL){ippsFree(fsect); fsect = NULL;}
}



bool MatPar::formCudaMatrixXY(){

	int volMatrix = na*nr;
	m_x = ippsMalloc_32f(volMatrix);
	m_y = ippsMalloc_32f(volMatrix);
	if(m_x == NULL || m_y == NULL){
		printf("Error: can not allocate memory (matrix).\n");
		return false;
	}

	Ipp32f *vcos, *vsin, *vq, *vecr;
	Ipp32f *vecx, *vecy;

	vcos = ippsMalloc_32f(na);
	vsin = ippsMalloc_32f(na);
	vq = ippsMalloc_32f(na);
	vecr = ippsMalloc_32f(nr);

	float dr, da;
	
	dr = (rmax-rmin)/float(nr-1);
	da = float(2.0*pi)/float(na-1);

	ippsVectorSlope_32f(vecr,nr,rmin,dr);
	ippsVectorSlope_32f(vq,na,0.0,da);
	ippsSinCos_32f_A24 (vq, vsin, vcos, na);
	ippsFree(vq); vq = NULL;

	vecx = ippsMalloc_32f(na);
	vecy = ippsMalloc_32f(na);

	float z, xc, ac, yc, ox, sdp;

	ac = 0.5f*float(ImageWidth) + AxisShift;
	xc = 0.5f*float(ImageWidth) + ConeCentreShiftX;
	yc = 0.5f*float(ImageHeight) + ConeCentreShiftY;
	ox = ac - xc;
	z = rowf-yc;
	sdp = SourceToDetector/PixelSize;
	//sdp = 40000.0/PixelSize;
	
	//printf("z: %f\n",z);
	//printf("rowf: %f\n",rowf);
	//printf("yc: %f\n",yc);
	for(int i=0;i<nr;i++){
		ippsMulC_32f(vcos,vecr[i],vecx,na);
		ippsMulC_32f(vsin,vecr[i],vecy,na);
		ippsAddC_32f_I(ox,vecx,na);
		ippsAddC_32f_I(sdp,vecy,na);
		ippsDivCRev_32f_I(sdp,vecy,na);
		ippsMul_32f(vecx,vecy,m_x+i*na,na);
		ippsMulC_32f(vecy,z,m_y+i*na,na);
	}
	ippsAddC_32f_I(xc,m_x,volMatrix);
	ippsAddC_32f_I(yc,m_y,volMatrix);
	
	ippsFree(vecx); vecx = NULL;
	ippsFree(vecy); vecy = NULL;
	ippsFree(vcos); vcos = NULL;
	ippsFree(vsin); vsin = NULL;
	ippsFree(vecr); vecr = NULL;

	ippsThreshold_LT_32f_I(m_x, volMatrix, 0.01f);
	ippsThreshold_GT_32f_I(m_x, volMatrix, float(ImageWidth) - 1.01f);

	ippsMinMax_32f(m_y,volMatrix,&y_min,&y_max);
	printf("m_y (before):\t%f\t%f\n",y_min,y_max);

	ippsThreshold_LT_32f_I(m_y, volMatrix, float(RestrMin)+1.01f);
	ippsThreshold_GT_32f_I(m_y, volMatrix, float(RestrMax)-1.01f);

	
	ippsMinMax_32f(m_y,volMatrix,&y_min,&y_max);
	printf("m_y:\t%f\t%f\n",y_min,y_max);
	y_mini = int(floor(y_min))-1;
	y_maxi = int(ceil(y_max))+1;
	
	printf("m_y (integer):\t%i\t%i\n",y_mini,y_maxi);

	ippsSubC_32f_I(float(y_mini),m_y,volMatrix);

	RowFirst = y_mini;
	

	return true;
}

