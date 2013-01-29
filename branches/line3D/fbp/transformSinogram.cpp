/* ----------------------------------------------------------------------
 * transformSinogram.cpp()
 * ----------------------------------------------------------------------
 * 
 * ----------------------------------------------------------------------
 */

#include "defs.h"

bool transform_sinogram_trans_fs(allData *aD){
	int nx, ny, mp, u, h, extra_left;
	int extra, tb_nx, tb_ny, ta_nx, ta_ny, orig_w;
	float da_orig, phi, ss, rad, a;
	float alpha, csc, sf, my1;
	xTransform *xt;
	xData *xd;
	xBackprojection *xb;

	nx = aD->nx;
	ny = aD->ny;
	alpha = aD->hms->fbp->inputData->pixelParam;
	xt = aD->hms->fbp->transform;
	xd = aD->data;
	xb = aD->hms->fbp->backprojection;

	if(!find_missed_projections(aD))return false;
	mp = aD->miss_proj;
	if(!check_crop(aD))return false;
	printTag(aD,"NumberOfPixelsToExtrapolateSinogram",xt->extrapolationPixels);
	if(xt->extrapolationPixels < 1 || xt->extrapolationPixels>aD->nx){
		printError(aD,"wrong number");
		return false;
	}
	if(xt->scaleType == SCALE_YES){
		printTag(aD,"ScaleWidth",xt->scaleWidth);
		if(xt->scaleWidth <1){
			printError(aD,"wrong number");
			return false;
		}
		printTag(aD,"ScaleHeight",xt->scaleHeight);
		if(xt->scaleHeight <1){
			printError(aD,"wrong number");
			return false;
		}
	}
	if(xt->extrapolationType != EXTRAPOLATION_NO){
		printTag(aD,"WidthAfterExtrapolation",(int)xt->extrapolationWidth);
		
		if(xt->scaleType == SCALE_YES){
			if(xt->extrapolationWidth < xt->scaleWidth){
				fprintf(aD->fp,"Error: wrong number (should be >= %u).\n",xt->scaleWidth);
				return false;
			}
		}else{
			if(xt->extrapolationWidth < aD->nx - xt->cropLeft - xt->cropRight){
				fprintf(aD->fp,"Error: wrong number (should be >= %u).\n",aD->nx);
				return false;
			}
		}
	}
		
	if(xt->rotationAngleType == ROTATION_ANGLE_OTHER){
		printTag(aD,"RotationAngle",xt->rotationAngle);
		if(xt->rotationAngle < 10.0 || xt->rotationAngle > 400.0){
			printError(aD,"wrong number");
			return false;
		}
	}
	tb_nx = nx+6;
	tb_ny = ny+aD->miss_proj+3;
	if(aD->hms->fbp->inputData->shape == SHAPE_PIXEL)tb_nx++;
	aD->tb_nx = tb_nx;
	aD->tb_ny = tb_ny;
	if(aD->hms->fbp->inputData->shape == SHAPE_PIXEL){
		xd->pp_a = ippsMalloc_64f(nx+1);
		xd->pp_b = ippsMalloc_64f(nx+1);
		xd->pp_c = ippsMalloc_64f(nx+1);
		xd->pp_f = ippsMalloc_64f(nx+1);
		xd->pp_t = ippsMalloc_64f(nx+1);
		if(xd->pp_a == NULL || xd->pp_b == NULL || xd->pp_c == NULL || xd->pp_f == NULL || xd->pp_t == NULL){
			fprintf(aD->fp,"Error: can not allocate memory (tridiagonal matrix).\n");
			return false;
		}
		if(alpha<1.e-6){
			fprintf(aD->fp,"Error: wrong value for the PixelParam: %f.\n",alpha);
			return false;
		}
		
	}

	if(xt->extrapolationType == EXTRAPOLATION_NO){
		if(xt->scaleType == SCALE_YES){
			ta_nx = xt->scaleWidth;
		}else{
			ta_nx = nx - xt->cropLeft - xt->cropRight;
			if(aD->hms->fbp->inputData->shape == SHAPE_PIXEL)ta_nx++;
		}
	}else{
		ta_nx = xt->extrapolationWidth;
	}

	if(xt->scaleType == SCALE_YES){
		ta_ny = xt->scaleHeight+1;
	}else{
		if(xt->rotationAngleType == ROTATION_ANGLE_180){
			ta_ny = ny + mp;
			if(xt->rotationAngleEndPoints == ROTATION_ANGLE_END_NO) ta_ny++;
		}else{
			h = ny + mp;
			if(xt->rotationAngleEndPoints == ROTATION_ANGLE_END_NO) h++;
			a = 180.0f*float(h)/xt->rotationAngle;
			ta_ny = int(a);
		}
	}

	xd->vtb = ippsMalloc_32f(aD->tb_nx*aD->tb_ny);
	if(check_is_null(aD, xd->vtb, "xd->vtb", "transform, before"))return false;
	
	aD->ta_nx = ta_nx;
	aD->ta_ny = ta_ny;
	xd->mapx = ippsMalloc_32f(ta_nx*ta_ny);
	xd->mapy = ippsMalloc_32f(ta_nx*ta_ny);
	xd->vta = ippsMalloc_32f(ta_nx*ta_ny);
	if(xd->mapx == NULL || xd->mapy == NULL || xd->vta == NULL){
		printError(aD,"can not allocate memory (transform, after)");
		return false;
	}

	ippsVectorSlope_32f(xd->mapx,ta_nx,0.0,1.0);
	for(int i=1;i<ta_ny;i++){
		ippsCopy_32f(xd->mapx,xd->mapx+i*ta_nx,ta_nx);
	}
	
		
	my1 = 0.0f;
	for(int i=0;i<ta_ny;i++){
		ippsSet_32f(float(i),xd->mapy+i*ta_nx,ta_nx);
	}
	
	u = ny+mp;
	if(xt->rotationAngleEndPoints == ROTATION_ANGLE_END_NO) u++;
	da_orig = float(u-1)/float(ta_ny-1);
	if(xt->rotationAngleType == ROTATION_ANGLE_OTHER)da_orig *= (180.0f/xt->rotationAngle);
	ippsMulC_32f_I(da_orig,xd->mapy,ta_nx*ta_ny);
	ippsAddC_32f_I(float(xt->cropBottom),xd->mapy,ta_nx*ta_ny);
	ippsThreshold_GTVal_32f_I(xd->mapy,ta_nx*ta_ny,float(ny+mp-1-xt->cropTop),float(tb_ny-1));

	extra = 0;
		
	orig_w = nx - xt->cropLeft - xt->cropRight;
	if(aD->hms->fbp->inputData->shape == SHAPE_PIXEL) orig_w++;
	if(xt->extrapolationType != EXTRAPOLATION_NO){
		if(xt->scaleType == SCALE_YES){
			extra = ta_nx - xt->scaleWidth;
		}else{
			extra = ta_nx - orig_w;
		}
	}
	extra_left = extra/2;
	ippsSubC_32f_I(float(extra_left),xd->mapx,ta_nx*ta_ny);
	
	csc = 1.0f;
		
	if(xt->scaleType == SCALE_YES){
		ippsMulC_32f_I(float(orig_w-1)/float(xt->scaleWidth-1),xd->mapx,ta_nx*ta_ny);
		csc = float(xt->scaleWidth-1)/float(orig_w-1);
	}
	aD->scale_factor = csc;

	//csh = float(-xt->cropLeft+extra_right);
	sf = extra_left+csc*(xb->imageCentre-xt->cropLeft);
	aD->ixc = (int)floor(sf);
	aD->dxc = (sf-aD->ixc)/csc;
	printTag(aD,"ImageCentreOriginal",xb->imageCentre);
	printTag(aD,"ImageCentreNew",aD->ixc);
		
	fprintf(aD->fp,"Re-center sinogram (radius): %f.\n",xt->reCentreRadius);
	if(fabs(xt->reCentreRadius)>float(aD->nx)){
		fprintf(aD->fp,"Error: wrong number.\n");
		return false;
	}
	if(fabs(xt->reCentreRadius)<1.e-3){
		aD->isReCentre = false;
	}else{
		aD->isReCentre = true;
	}
	fprintf(aD->fp,"Re-center sinogram (angle): %f.\n",xt->reCentreAngle);
	xt->reCentreAngle *= float(pi/180.0);

	ss = acos(-1.0f)/float(ta_ny-1);
	for(int i=0;i<ta_ny;i++){
		phi = xt->reCentreAngle + float(i)*ss;
		rad = xt->reCentreRadius * cos(phi);
		ippsSubC_32f_I(rad,xd->mapx+i*ta_nx,ta_nx);
	}
	ippsAddC_32f_I(float(xt->cropLeft+3+aD->dxc),xd->mapx,ta_nx*ta_ny);
	if(aD->hms->fbp->outputData->type == OUTPUT_TYPE_SOLUTION){
		for(int i=0;i<ta_ny;i++){
			sf = float(i) + my1;
			if(!(sf<0.0 || sf>float(ta_ny-1)))continue;
			ippsSubCRev_32f_I(xb->imageCentre,xd->mapx+i*ta_nx,ta_nx);
			ippsAddC_32f_I(xb->imageCentre,xd->mapx+i*ta_nx,ta_nx);
		}
	}
		
	ippsThreshold_LTVal_32f_I(xd->mapx,ta_nx*ta_ny,float(xt->cropLeft+3),1.0f);
	ippsThreshold_GTVal_32f_I(xd->mapx,ta_nx*ta_ny,float(nx+2-xt->cropRight),float(tb_nx)-2.0f);
	return true;
}




bool transform_sinogram_sol_fs(allData *aD){
	char *block_wh;
	unsigned int nx, ny, mp, h, extra_left, extra_right, extra;
	unsigned int tb_nx, tb_ny, ta_nx, ta_ny;
	unsigned int ixc, rr;
	unsigned int dwo, dho, optw, opth, optwh;
	unsigned int diff, rval, nz;
	
	float a, alpha;
	float new_xc, shift, fxc, dxc, xc, da;
	
	xTransform *xt;
	xData *xd;
	xRoi *xr;
	xBackprojection *xb;
	gpu_info *gi;

	block_wh = NULL;

	nx = aD->nx;
	ny = aD->ny;
	alpha = aD->hms->fbp->inputData->pixelParam;
	xt = aD->hms->fbp->transform;
	xd = aD->data;
	xr = aD->hms->fbp->backprojection->roi;
	xb = aD->hms->fbp->backprojection;

	gi = aD->gi;


	if(!get_gpu_info(aD))return false;
	if(!check_crop(aD))return false;
	if(!check_roi(aD))return false;

	gi->xmin = xr->xmin;
	gi->xmax = xr->xmax;
	gi->ymin = xr->ymin;
	gi->ymax = xr->ymax;

	if(!check_roi_size(aD))return false;

	if(!find_missed_projections(aD))return false;
	mp = aD->miss_proj;
	
	h = ny + mp;
	if(xt->rotationAngleEndPoints == ROTATION_ANGLE_END_NO) h++;

	if(xt->rotationAngleType == ROTATION_ANGLE_180){
		ta_ny = h;
	}else{
		a = 180.0f*float(h)/xt->rotationAngle;
		ta_ny = int(a);
	}

	da = pif/float(ta_ny-1);
	gi->rotAngleStep = da;
	
	optwh = gi->maxResidentThreads / gi->maxResidentBlocks;
	printTag(aD,"OptimalNumberOfThreadsInBlock",optwh);
	diff = optwh;
	optw = optwh;
	opth = 1;
	for(unsigned int z = 8; z<optwh; z+= 8){
		if(optwh/z == 1)continue;
		rval = (unsigned int)abs((int)(z-(optwh/z)));
		if(rval > diff) continue;
		optw = z;
		opth = optwh/z;
		diff = rval;
	}

	printTag(aD,"OptimalGPUBlockWidth",optw);
	printTag(aD,"OptimalGPUBlockHeight",opth);

	dwo = (int)(DEFAULT_GPU_BLOCK_WIDTH);
	dho = (int)(DEFAULT_GPU_BLOCK_HEIGHT);
	printTag(aD,"DefaultGPUBlockWidth",dwo);
	printTag(aD,"DefaultGPUBlockHeight",dho);

	block_wh = getenv("GPU_BLOCK_WIDTH");
	if(block_wh == NULL){
		printInfo(aD,"No GPU_BLOCK_WIDTH is given");
		dwo = DEFAULT_GPU_BLOCK_WIDTH;
	}else{
		dwo = atoi(block_wh);
		printTag(aD,"GpuBlockWidthSys",dwo);
	}

	block_wh = getenv("GPU_BLOCK_HEIGHT");
	if(block_wh == NULL){
		printInfo(aD,"No GPU_BLOCK_HEIGHT is given");
		dho = DEFAULT_GPU_BLOCK_HEIGHT;
	}else{
		dho = atoi(block_wh);
		printTag(aD,"GpuBlockHeightSys",dwo);
	}

	if(dwo <= 0 || dho <= 0 || dwo * dho > gi->maxThreadsPerBlock){
		dwo = optw;
		dho = opth;
	}

	//dwo = 8;
	//dho = 8;

	printTag(aD,"GPUBlockWidth",dwo);
	printTag(aD,"GPUBlockHeight",dho);

	gi->wo = dwo;
	gi->ho = dho;
	find_params(aD);
	ta_nx = gi->flength;

	printTag(aD,"NumberOfPixelsToExtrapolateSinogram",(int)xt->extrapolationPixels);
	if(xt->extrapolationPixels < 1 || xt->extrapolationPixels > aD->nx){
		printError(aD,"wrong number");
		return false;
	}
	printInfo(aD,"No scaling");
	
	if(xt->rotationAngleType == ROTATION_ANGLE_OTHER){
		printTag(aD,"RotationAngle",xt->rotationAngle);
		if(xt->rotationAngle < 10.0 || xt->rotationAngle > 400.0){
			printError(aD,"wrong number");
			return false;
		}
	}
	tb_nx = nx+6;
	tb_ny = ny+aD->miss_proj+3;
	if(aD->hms->fbp->inputData->shape == SHAPE_PIXEL)tb_nx++;
	aD->tb_nx = tb_nx;
	aD->tb_ny = tb_ny;
	if(aD->hms->fbp->inputData->shape == SHAPE_PIXEL){
		xd->pp_a = ippsMalloc_64f(nx+1);
		xd->pp_b = ippsMalloc_64f(nx+1);
		xd->pp_c = ippsMalloc_64f(nx+1);
		xd->pp_f = ippsMalloc_64f(nx+1);
		xd->pp_t = ippsMalloc_64f(nx+1);
		if(xd->pp_a == NULL || xd->pp_b == NULL || xd->pp_c == NULL || xd->pp_f == NULL || xd->pp_t == NULL){
			printError(aD,"can not allocate memory (tridiagonal matrix)");
			return false;
		}
		if(alpha<1.e-6){
			sprintf(aD->message,"wrong value for PixelParam: %f",alpha);
			printError(aD);
			return false;
		}
	}
	aD->ta_ny_13 = ta_ny;
	int_ratio(ta_ny, gi->vertH, &rr);
	ta_ny = rr;

	gi->nxo = ta_nx;
	gi->nyo = ta_ny;

	xd->vtb = ippsMalloc_32f(aD->tb_nx*aD->tb_ny);
	if(check_is_null(aD, xd->vtb, "xd->vtb", "transform, before"))return false;

	aD->nz = 4;	
	nz = aD->nz;

	aD->ta_nx = ta_nx;
	aD->ta_ny = ta_ny;
	xd->mapx = ippsMalloc_32f(ta_nx*ta_ny);
	if(check_is_null(aD, xd->mapx, "xd->mapx", "transform, after"))return false;
	xd->mapy = ippsMalloc_32f(ta_nx*ta_ny);
	if(check_is_null(aD, xd->mapy, "xd->mapy", "transform, after"))return false;
	xd->vta = ippsMalloc_32f(ta_nx*ta_ny * nz);
	if(check_is_null(aD, xd->vta, "xd->vta", "transform, after"))return false;
	

	ippsVectorSlope_32f(xd->mapx,ta_nx,0.0,1.0);
	for(unsigned int i=1;i<ta_ny;i++){
		ippsCopy_32f(xd->mapx,xd->mapx+i*ta_nx,ta_nx);
	}

	extra = gi->flength - gi->origw;
	extra_left = extra/2;
	extra_right = extra - extra_left;

	xc = xb->imageCentre;
	fxc = float(extra_left) + xc - float(xt->cropLeft);
	ixc = (int)(ceil(fxc/float(gi->warpSize)))*gi->warpSize;
	dxc = float(ixc)-fxc;
	if(dxc > float(extra_right))dxc = 0.0f;

	
	new_xc = float(extra_left) + xc + dxc;
	dxc = new_xc - xc;
	aD->new_xc = new_xc;

	printTag(aD,"ImageCentreOriginal",xc);
	printTag(aD,"ImageCentreNew",new_xc);

	
	for(unsigned int i=0;i<ta_ny;i++){
		shift = dxc - gi->roiR*cosf(gi->roiA-da*float(i)) - 3.0f;
		ippsSubC_32f_I(shift,xd->mapx+i*ta_nx,ta_nx);
	}	
	
	for(unsigned int i=0;i<ta_ny;i++){
		ippsSet_32f(float(i),xd->mapy+i*ta_nx,ta_nx);
	}

	ippsAddC_32f_I(float(xt->cropBottom),xd->mapy,ta_nx*ta_ny);
	ippsThreshold_GTVal_32f_I(xd->mapy,ta_nx*ta_ny,float(ny+mp-1-xt->cropTop),float(tb_ny-1));
			
	ippsThreshold_LTVal_32f_I(xd->mapx,ta_nx*ta_ny,float(xt->cropLeft+3),1.0f);
	ippsThreshold_GTVal_32f_I(xd->mapx,ta_nx*ta_ny,float(nx+2-xt->cropRight),float(tb_nx)-2.0f);
	return true;
}


bool transform_sinogram(allData *aD){
	int nx, ny, mp;
	int tb_nx, tb_ny, tb_nt;
	int ta_nx, ta_ny;
	int i1, i2, i3;
	int k;
	int s;
	float w;
	Ipp32f sv;
	double mv;
	float alpha;
	xTransform *xt;
	xData *xd;

	nx = aD->nx;
	ny = aD->ny;
	alpha = aD->hms->fbp->inputData->pixelParam;
	xt = aD->hms->fbp->transform;
	xd = aD->data;

	printTagStart(aD,"Transform","Transform the sinogram");

	if(aD->isFirstSlice){
		if(aD->hms->fbp->outputData->type == OUTPUT_TYPE_SOLUTION){
			printInfo(aD,"Here1");
			if(!transform_sinogram_sol_fs(aD))return false;
			printInfo(aD,"Here2");
		}else{
			if(!transform_sinogram_trans_fs(aD))return false;
		}
	}

	printInfo(aD,"Here0");
	
	mp = aD->miss_proj;
	
	tb_nx = aD->tb_nx;
	tb_ny = aD->tb_ny;
	tb_nt = tb_nx*tb_ny;
	ta_nx = aD->ta_nx;
	ta_ny = aD->ta_ny;

	ippsZero_32f(xd->vtb,tb_nt);
	
	i1 = -1;
	k = 0;
	for(int i=0;i<mp+1;i++){
		i2 = xd->vmiss[i];
		if(i2 - i1 == 1){
			k++;
			i1 = i2;
			continue;
		}
		for(int j=i1-k+1;j<i2-k;j++){
			
			if(aD->hms->fbp->inputData->shape == SHAPE_PIXEL){
				ippsConvert_32f64f(xd->veci+j*nx,xd->pp_t,nx);
				ippsSet_64f(2.0+alpha,xd->pp_b,nx);
				ippsSet_64f(1.0-alpha,xd->pp_a,nx+1);
				ippsSet_64f(1.0-alpha,xd->pp_c,nx+1);
				ippsZero_64f(xd->pp_f,nx+1);
				ippsAdd_64f_I(xd->pp_t,xd->pp_f,nx);
				ippsAdd_64f_I(xd->pp_t,xd->pp_f+1,nx);
				xd->pp_b[0] = 1.0+alpha;
				//xd->pp_b[1] = 2.0+alpha;
				//xd->pp_a[1] = 1.0-alpha;
				//xd->pp_c[0] = 1.0-alpha;

				xd->pp_b[nx] = 1.0+alpha;
				//xd->pp_b[nx-1] = 2.0+alpha;
				//xd->pp_a[nx] = 1.0-alpha;
				//xd->pp_c[nx-1] = 1.0-alpha;

				xd->pp_c[0]/=xd->pp_b[0];
				xd->pp_f[0]/=xd->pp_b[0];
				for (int jj = 1; jj <= nx; jj++){
					mv = xd->pp_b[jj] - xd->pp_c[jj-1] * xd->pp_a[jj];
					xd->pp_c[jj]/=mv;
					xd->pp_f[jj]=(xd->pp_f[jj] - xd->pp_f[jj-1]*xd->pp_a[jj])/mv;
				}

				xd->pp_t[nx] = xd->pp_f[nx];
			 
				for (int jj = nx - 1; jj >= 0; jj--){
					xd->pp_t[jj]=xd->pp_f[jj]-xd->pp_c[jj]*xd->pp_t[jj+1];
				}

				ippsConvert_64f32f(xd->pp_t,xd->vtb+(j+k)*tb_nx+3,nx+1);
			}else{
				ippsCopy_32f(xd->veci+j*nx,xd->vtb+(j+k)*tb_nx+3,nx);
			}
		}
		k++;
		i1=i2;
	}

	if(xt->missedProjectionsType == MISSED_PROJECTIONS_LINEAR && mp >0){
		for(int i = 0;i<mp;i++){
			i3 = xd->vmiss[i];
			i1 = i3-1;
			i2 = i3+1;
			s = 1;
			for(;;){
				if(i1<0)break;
				if(xd->vmiss[i-s] != i1)break;
				i1--;
				s++;
			}
			s = 1;
			for(;;){
				if(i2 >= (int)aD->ny)break;
				if(xd->vmiss[i+s] != i2)break;
				i2++;
				s++;
			}
			if(i1 <0) i1 = i2;
			if(i2>= (int)aD->ny) i2 = i1;

			//printf("DDD: %i, %i\n",i1,i2);
			if(i1 == i2){
				ippsCopy_32f(xd->vtb+i1*tb_nx,xd->vtb+i3*tb_nx,tb_nx);
			}else{
				w = float(i3-i1)/float(i2-i1);
				ippsAddProductC_32f(xd->vtb+i1*tb_nx,1.0f-w,xd->vtb+i3*tb_nx,tb_nx);
				ippsAddProductC_32f(xd->vtb+i2*tb_nx,w,xd->vtb+i3*tb_nx,tb_nx);

			}
		}
	}
	
	for(int i=0;i<ny+mp;i++){
		ippsMean_32f(xd->vtb+i*tb_nx+3+xt->cropLeft,xt->extrapolationPixels,&sv,ippAlgHintNone);
		ippsSet_32f(sv,xd->vtb+i*tb_nx,3);
		ippsMean_32f(xd->vtb+i*tb_nx+aD->nx-xt->extrapolationPixels+3-xt->cropRight,xt->extrapolationPixels,&sv,ippAlgHintNone);
		ippsSet_32f(sv,xd->vtb+(i+1)*tb_nx-3,3);
	}
	IppiRect srcRoi = {0,0,tb_nx,tb_ny};
	IppiSize srcSize = {tb_nx,tb_ny};
	IppiSize dstRoiSize = {ta_nx,ta_ny};
	
	ippiRemap_32f_C1R(xd->vtb, srcSize, 4*tb_nx, srcRoi, xd->mapx, 4*ta_nx, xd->mapy, 4*ta_nx, xd->vta, 4*ta_nx, dstRoiSize, xt->interpolation);

	int nz;
	nz = aD->nz;
	for(int i=1;i<nz;i++){
		ippsCopy_32f(xd->vta,xd->vta+i*ta_nx*ta_ny,ta_nx*ta_ny);
	}

	printInfo(aD,"Here");
	printTagEnd(aD);

	return true;
}



