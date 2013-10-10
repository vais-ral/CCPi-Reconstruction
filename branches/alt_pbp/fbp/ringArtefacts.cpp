/* ----------------------------------------------------------------------
 * ringArtefacts.cpp()
 * ----------------------------------------------------------------------
 * Suppression of ring artefacts
 * ----------------------------------------------------------------------
 */

#include "defs.h"

bool ring_artefacts_removal_fs(allData *aD){
	
	int nx, nxi, ny, ii, jj;
	int t, cL, cR, nS;
	Ipp64f mdiv, alpha, tau, paramN, paramR;
	Ipp64f *vec_exp, *vec1, *vec2;
	xData *xd;

	nxi = aD->nx;
	ny = aD->ny;
	xd = aD->data;

	t = aD->hms->fbp->preprocessing->ringArtefacts->type;
	
	if(!check_crop(aD))return false;
	
	cL = aD->hms->fbp->transform->cropLeft;
	cR = aD->hms->fbp->transform->cropRight;

	nx = nxi - (cL + cR);
			
	paramN = double(aD->hms->fbp->preprocessing->ringArtefacts->parameterN);
	paramR = double(aD->hms->fbp->preprocessing->ringArtefacts->parameterR);
	

	printTag(aD,"Parameter_N",float(paramN),"ratio of norms");
	printTag(aD,"Parameter_R",float(paramR),"regularisation");
	
	if(paramN <-1e-10 || paramR <1e-8){
		printError(aD,"wrong parameters");
		return false;
	}

	nS = aD->hms->fbp->preprocessing->ringArtefacts->num_series;
	printTag(aD,"NumSeries",nS,"number of series");
	if(nS<1 || nS >100){
		printError(aD,"wrong parameter");
		return false;
	}
	
	if(t == RING_ARTEFACTS_COLUMN || t == RING_ARTEFACTS_AML){
		xd->vec32 = ippsMalloc_32f(nx);
		xd->vec64 = ippsMalloc_64f(nx);
		xd->vec_res = ippsMalloc_64f(nx);
		if(xd->vec32 == NULL || xd->vec64 == NULL || xd->vec_res == NULL){
			printError(aD,"can not allocate memory (ring artefacts)");
			return false;
		}
	}

	

	if(t == RING_ARTEFACTS_AML){
		xd->vecRS = ippsMalloc_64f(2*nS*ny);
		vec1 = ippsMalloc_64f(1);
		vec2 = ippsMalloc_64f(1);
		vec_exp = ippsMalloc_64f(2*nx+2);
		xd->mata = ippsMalloc_64f(nS*nx*nx);
		if(vec1 == NULL || vec2 == NULL || xd->mata == NULL || vec_exp == NULL || xd->vecRS == NULL){
			printError(aD,"can not allocate memory (ring artefacts, matrix)");
			return false;
		}

		
		for(int s = 1; s <nS; s++){
			ippsVectorSlope_64f(xd->vecRS,ny,0.0,double(s)*pi/float(ny));
			ippsSinCos_64f_A53(xd->vecRS,xd->vecRS+(2*s-1)*ny,xd->vecRS+2*s*ny,ny);
			ippsMulC_64f_I(sqrt(2.0),xd->vecRS+(2*s-1)*ny,2*ny);
		}
		ippsSet_64f(1.0,xd->vecRS,ny);
	
		for(int s = 1; s<=nS; s++){
			alpha = (paramN+paramR)*double(s*s);
			mdiv = sqrt(alpha*(alpha+4.0));
			vec1[0]=sqrt(alpha)/2.0;
			ippsAsinh_64f_A53 (vec1, vec2, 1);
			tau = 2.0*vec2[0];
			ippsVectorSlope_64f(vec_exp,2*nx+2,0.0,-tau);

			ippsExp_64f_I(vec_exp, 2*nx+2);

			for(int i = 0;i<nx;i++){
				for(int j=0;j<nx;j++){
					if(j>i){
						ii=i;
						jj=j;
					}else{
						ii=j;
						jj=i;
					}
					ii++;
					jj++;
					xd->mata[i*nx+j+(s-1)*nx*nx] = vec_exp[abs(jj-ii)]* (1.0+vec_exp[2*nx-2*jj+1]) * (1.0+vec_exp[2*ii-1]);
				}
			}
			ippsDivC_64f_I(mdiv*(1.0-vec_exp[2*nx]),xd->mata+(s-1)*nx*nx,nx*nx);
		}

		ippsFree(vec1); vec1 = NULL;
		ippsFree(vec2); vec2 = NULL;
		ippsFree(vec_exp); vec_exp = NULL;
			
	}
	

	
	return true;
}


bool ring_artefacts_removal(allData *aD){
	printTagStart(aD,"RingArtefacts");

	int nx, nxi, ny;
	int t, cL, cR, nS;
	Ipp64f paramN;
	xData *xd;

	nxi = aD->nx;
	ny = aD->ny;
	xd = aD->data;

	t = aD->hms->fbp->preprocessing->ringArtefacts->type;
	
	paramN = double(aD->hms->fbp->preprocessing->ringArtefacts->parameterN);
//	paramR = double(aD->hms->fbp->preprocessing->ringArtefacts->parameterR);

	nS = aD->hms->fbp->preprocessing->ringArtefacts->num_series;

	cL = aD->hms->fbp->transform->cropLeft;
	cR = aD->hms->fbp->transform->cropRight;

	nx = nxi - (cL + cR);
	
	if(aD->isFirstSlice){
		if(!ring_artefacts_removal_fs(aD))return false;
	}

	int ss;

/*	ippsZero_32f(xd->veci,nxi*ny);

	for(int i=ny/3;i<(0.6*ny);i++){
		xd->veci[i*nxi+2000]+= 100.0;//(10.0*cos(pi*i/float(ny)));

	//	xd->veci[i*nxi+1500]+= (10.0*sin(10.5*pi*i/float(ny)));
	}
	
*/	switch(t){
		case RING_ARTEFACTS_NO:
			printInfo(aD,"No ring artefact suppression");
			break;
		case RING_ARTEFACTS_COLUMN:
			printInfo(aD,"A mean value is subtracted from each column");
			ippsZero_64f(xd->vec_res,nx);
			for(int i=0;i<ny;i++){
				ippsConvert_32f64f(xd->veci+i*nxi + cL,xd->vec64,nx);
				ippsAdd_64f_I(xd->vec64,xd->vec_res,nx);
			}
			ippsDivC_64f_I(double(ny),xd->vec_res,nx);
			ippsConvert_64f32f(xd->vec_res,xd->vec32,nx);
			for(int i=0;i<ny;i++){
				ippsSub_32f_I(xd->vec32,xd->veci + i*nxi + cL,nx);
			}
			break;
		case RING_ARTEFACTS_AML:
			printInfo(aD,"Suppression based on paper \"Applied Mathematics Letters, v. 23, p. 1489\"");

			for(int s = 0; s < 2*nS-1; s++){
				ss = (s+2)/2;
				
				ippsZero_64f(xd->vec_res,nx);
					
				for(int i=0;i<ny;i++){
					ippsConvert_32f64f(xd->veci+i*nxi + cL,xd->vec64,nx);
					ippsMulC_64f_I(xd->vecRS[s*ny+i],xd->vec64,nx);
					ippsAdd_64f_I(xd->vec64,xd->vec_res,nx);
				}
				ippsDivC_64f_I(double(ny),xd->vec_res,nx);
			

				ippsZero_64f(xd->vec64,nx);
				//ippsMulC_64f(xd->vec_res,paramN*paramR*double(ss*ss),xd->vec64,nx);
				ippsMulC_64f(xd->vec_res,paramN*double(ss*ss),xd->vec64,nx);
				ippsAdd_64f_I(xd->vec_res,xd->vec64,nx-1);
				ippsSub_64f_I(xd->vec_res+1,xd->vec64,nx-1);
				ippsAdd_64f_I(xd->vec_res+1,xd->vec64+1,nx-1);
				ippsSub_64f_I(xd->vec_res,xd->vec64+1,nx-1);

				for(int i=0;i<nx;i++){
					ippsDotProd_64f(xd->vec64,xd->mata+i*nx,nx,xd->vec_res + i);
				}
				ippsConvert_64f32f(xd->vec_res,xd->vec32,nx);
			
				for(int i=0;i<ny;i++){
					//ippsSub_32f_I(xd->vec32,xd->veci+i*nxi + cL,nx);
					ippsAddProductC_32f(xd->vec32,-float(xd->vecRS[s*ny+i]),xd->veci+i*nxi + cL,nx);
				}


			/*	ippsZero_64f(xd->vec_res,nx);
					
				for(int i=0;i<ny;i++){
					ippsConvert_32f64f(xd->veci+i*nxi + cL,xd->vec64,nx);
					ippsMulC_64f_I(xd->vecRS[s*ny+i],xd->vec64,nx);
					ippsAdd_64f_I(xd->vec64,xd->vec_res,nx);
				}
				ippsDivC_64f_I(double(ny),xd->vec_res,nx);
				ss = 0.0;*/
			}
			break;
		default:
			break;
	}

	printTagEnd(aD);
	return true;
}


