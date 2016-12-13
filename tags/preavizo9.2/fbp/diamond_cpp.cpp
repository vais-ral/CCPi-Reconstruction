#include "defs.h"



bool form_filter(allData *aD){
	int nx;
	int n1, n2;
	int ns, nl;
	float omega, pixel_size;
	float fxc, frad;
	float gamma, sg, x, su;
	Ipp32f sr, so;
	Ipp32f *v1, *v2, *v3, *v4;
	xData *xd;
	xFilter *xf;

	printTagStart(aD,"FilterForFBP","Forming a filter for FBP");

	xf = aD->hms->fbp->backprojection->filter;
	pixel_size = xf->pixelSize;
	printTag(aD,"PixelSizeInput",pixel_size);
	
	if(pixel_size <1.e-6){
		printError(aD,"wrong number");
		return false;
	}
	
	printTag(aD,"PixelSizeOutput",pixel_size);
	printInfo(aD);

	xd = aD->data;
	
	nx = aD->ta_nx;
	xd->vfilter = ippsMalloc_32f(nx);
	
	v1 = ippsMalloc_32f(nx);
	v2 = ippsMalloc_32f(nx);
	v3 = ippsMalloc_32f(nx);
	v4 = ippsMalloc_32f(nx);
	if(v1 == NULL || v2 == NULL || v3 == NULL || v4 == NULL || xd->vfilter == NULL){
		printError(aD,"can not allocate memory (filter)");
		return false;
	}

	omega = xf->bandwidth;
	printTag(aD,"FilterBandwidth",omega);
	if(omega > 1.0f || omega < 0.00001f){
		printError(aD,"wrong number (should be in (0, 1])");
		return false;
	} 

	ippsVectorSlope_32f(v1,nx,0,1.0);
	ippsThreshold_GTVal_32f(v1,v4,nx,float(nx)*omega,0.0);

	switch(xf->name){
		case FILTER_NAME_RL:
			printTag(aD,"FilterName", "Ram-Lak");
			ippsCopy_32f(v4,v1,nx);
			
			break;
		case FILTER_NAME_SL:
			printTag(aD,"FilterName", "Shepp-Logan");
			ippsMulC_32f_I(pif/(float(nx)*omega),v1,nx);
			ippsThreshold_GTVal_32f_I(v1,nx,0.5f*pif,0.0f);
			for(int i=0;i<nx;i++){
				v2[i] = sinc(v1[i]);
			}
			ippsMul_32f(v2,v4,v1,nx);
			break;
		case FILTER_NAME_COS:
			printTag(aD,"FilterName", "Cosine");
			ippsMulC_32f_I(pif/(float(nx)*omega),v1,nx);
			ippsThreshold_GTVal_32f_I(v1,nx,0.5f*pif,0.0f);
			ippsCos_32f_A24(v1,v2,nx);
			ippsMul_32f(v2,v4,v1,nx);
			break;
		default:
			printError(aD,"unknown filter");
			return false;
			break;
	}
	ippsCopy_32f(v1,v4,nx);

	ippsVectorSlope_32f(v1,nx,0,1.0f);
	ippsMulC_32f_I(2.0f*pif/float(nx),v1,nx);
	ippsMulC_32f(v1,2.0,v2,nx);
	ippsCos_32f_A24(v1,v3,nx);
	ippsCopy_32f(v3,v1,nx);
	ippsCos_32f_A24(v2,v3,nx);
	 
	switch(xf->windowName){
		case FILTER_WINDOW_NO:
			printTag(aD,"Window", "No");
			ippsSet_32f(1.0,v3,nx);
			break;
		case FILTER_WINDOW_HANN:
			printTag(aD,"Window", "Hann");
			ippsMulC_32f_I(0.50f,v1,nx);
			ippsAddC_32f(v1,0.50f,v3,nx);
			break;
		case FILTER_WINDOW_HAMMING:
			printTag(aD,"Window", "Hamming");
			ippsMulC_32f_I(0.46f,v1,nx);
			ippsAddC_32f(v1,0.54f,v3,nx);
			break;
		case FILTER_WINDOW_BLACKMAN:
			printTag(aD,"Window", "Blackman");
			ippsMulC_32f_I(0.08f,v3,nx);
			ippsAddC_32f_I(0.42f,v3,nx);
			ippsAddProductC_32f(v1,0.5,v3,nx);
			break;
		default:
			printError(aD,"unknown window");
			return false;
			break;
	}
	ippsMul_32f_I(v3,v4,nx);

	
	n1 = (nx+1)/2;
	n2 = nx - n1;
	ippsCopy_32f(v4,xd->vfilter,n1);
	ippsFlip_32f(v4+1,xd->vfilter+n1,n2);

	ippsFree(v1); v1 = NULL;
	ippsFree(v2); v2 = NULL;
	ippsFree(v3); v3 = NULL;
	ippsFree(v4); v4 = NULL;

	
	xd->vecti = ippsMalloc_32f(nx);
	xd->vecto = ippsMalloc_32f(nx);
	xd->vectr = ippsMalloc_32f(nx);
	if(xd->vecti == NULL || xd->vecto == NULL || xd->vectr == NULL){
		printError(aD,"can not allocate memory (normalisation)");
		return false;
	}

	
	fxc = float(nx-1)*0.5f;
	frad = 0.1f*float(nx-1);
	ippsVectorSlope_32f(xd->vecti, nx, -fxc, 1.0f);
	ippsAbs_32f_I(xd->vecti, nx);

	ippsCopy_32f(xd->vecti, xd->vectr, nx);
	ippsCopy_32f(xd->vecti, xd->vecto, nx);

	ippsSqr_32f_I(xd->vecti, nx);
	ippsSubCRev_32f_I(frad*frad,xd->vecti,nx);
	ippsThreshold_LT_32f_I(xd->vecti,nx,0.0f);
	ippsSqrt_32f_I(xd->vecti,nx);
	
	

	switch(xf->norm){
		case FILTER_NORM_CC: 
		case FILTER_NORM_CA:
			ippsThreshold_LTValGTVal_32f_I(xd->vectr,nx,frad,1.0,frad,0.0);
			ippsMulC_32f_I(2.0,xd->vecti,nx);
			break;
		case FILTER_NORM_PC: 
		case FILTER_NORM_PA:
			ippsSubCRev_32f_I(frad,xd->vectr,nx);
			ippsThreshold_LT_32f_I(xd->vectr,nx,0.0);
			for(int i=0;i<nx;i++){
				x = xd->vecto[i];
				gamma = x/frad;
				sg = 1-gamma*gamma;
				if(sg<=0.0){
					xd->vecti[i]=0.0;
					continue;
				}
				sg = sqrt(sg);
				xd->vecti[i] = frad*frad*(sg-gamma*gamma*log((1.0f+sg)/gamma));
			}
			//ippsSubCRev_32f_I(frad,xd->vectr,nx);
			//ippsThreshold_LT_32f_I(xd->vectr,nx,0.0);
			//ippsMul_32f_I(xd->vectr,xd->vecti,nx);
			break;
		case FILTER_NORM_SC: 
		case FILTER_NORM_SA:
			ippsCopy_32f(xd->vecti,xd->vectr,nx);
			ippsSqr_32f_I(xd->vecti,nx);
			ippsMulC_32f_I(0.5f*pif,xd->vecti,nx);
			break;
		default:
			printError(aD,"unknown sample for normalisation");
			return false;
			break;
	}

	Ipp32f *zero = ippsMalloc_32f(nx);
	xd->veccF = ippsMalloc_32fc(nx);
	xd->veccI = ippsMalloc_32fc(nx);
	ippsZero_32f(zero, nx);
	ippsRealToCplx_32f(xd->vfilter, zero, xd->veccF, nx);
	ippsRealToCplx_32f(xd->vecti, zero, xd->veccI, nx);
	
	ippsFree(zero); zero = NULL;

	if(!fbp_axial_cuda(aD))return false;

	ippsFree(xd->veccF); xd->veccF = NULL;
	ippsFree(xd->veccI); xd->veccI = NULL;

	ns = int(0.47f*float(nx));
	nl = int(0.06f*float(nx));

	sr = 1.0f;
	so = 1.0f;
	
	switch(xf->norm){
		case FILTER_NORM_CC:
		case FILTER_NORM_PC: 
		case FILTER_NORM_SC: 
			ippsMean_32f(xd->vectr+ns,nl,&sr,ippAlgHintNone);
			ippsMean_32f(xd->vecto+ns,nl,&so,ippAlgHintNone);
			break;
		case FILTER_NORM_CA:
		case FILTER_NORM_PA: 
		case FILTER_NORM_SA: 
			ippsMean_32f(xd->vectr,nx,&sr,ippAlgHintNone);
			ippsMean_32f(xd->vecto,nx,&so,ippAlgHintNone);
			break;
		default:
			break;
	}

	su = sr/(so*pixel_size);
	ippsMulC_32f_I(su,xd->vfilter,nx);

	ippsFree(xd->vecti); xd->vecti = NULL;
	ippsFree(xd->vecto); xd->vecto = NULL;
	ippsFree(xd->vectr); xd->vectr = NULL;

	printTagEnd(aD);
	return true;
}

bool check_rows(allData *aD){
	aD->row_first = aD->hms->fbp->inputData->fileFirst;
	aD->row_last = aD->hms->fbp->inputData->fileLast;
	aD->row_step = aD->hms->fbp->inputData->fileStep;

	if(aD->row_step <1){
		sprintf(aD->message,"the step for slices (%u, should be >0)",aD->row_step);
		printError(aD);
		return false;
	}

	if(aD->row_last < aD->row_first){
		sprintf(aD->message,"the first slice (%u) > the last slice (%u)",aD->row_first, aD->row_last);
		printError(aD);
		return false;
	}

	sprintf(aD->message,"The slices: first (%u), last (%u), step (%u)",aD->row_first,aD->row_last,aD->row_step);
	printInfo(aD);

	unsigned int r;
	for(r = aD->row_first; r<= aD->row_last; r+= aD->row_step){}
	r-=aD->row_step;
	if(aD->row_last != r){
		aD->row_last = r;
		sprintf(aD->message,"The last slice will be %u",r);
		printInfo(aD);
	}

	if(aD->hms->fbp->flatDarkFields->darkField->type == FIELD_TYPE_USER && aD->hms->fbp->flatDarkFields->flatField->type == FIELD_TYPE_USER) return true;

	if(aD->row_last >= aD->dark_ny){
		sprintf(aD->message,"the last slice (%u) >= the last row of dark/flat field (%u)",aD->row_last, aD->dark_ny);
		printError(aD);
		return false;
	}

	if(aD->row_last >= aD->flat_ny){
		sprintf(aD->message,"the last slice (%u) >= the last row of dark/flat field (%u)",aD->row_last, aD->flat_ny);
		printError(aD);
		return false;
	}

	return true;
}




allData::allData(xHMset *hmset){
	cspace[0] = '\0';
	fp = NULL;
	tagNumber = 0;
	gi = new gpu_info();
	data = new xData();
	time = new xTime();
	fg = new xFBP_gpu();
	hms = hmset;
	isFirstSlice = true;
	//dataSource = NULL;
	//vecDark = NULL;
	//vecLight = NULL;
	//matrDark1 = NULL;
	//matrDark2 = NULL;
}

allData::~allData(){
	delete fg; fg = NULL;
	delete gi; gi = NULL;
	delete data; data = NULL;
	delete time; time = NULL;
	hms = NULL;
}


bool allocate_memory(allData *aD){
	float fm;
	int m, mm;
	xData *xd;
	xInputData *xid;

	xd = aD->data;
	xid = aD->hms->fbp->inputData;

	m = xid->memorySizeMin;
	printTag(aD,"MinimalSizeOfFiles",m,"in bytes");
	if(m<0){
		printError(aD,"negative value");
		return false;
	}
	fm = xid->memorySizeFMax;
	mm = (int)(1024.0f*1024.0f*fm);

	printTag(aD,"MaximalSizeOfFiles",fm,"in MB");
	printTag(aD,"MaximalSizeOfFiles",mm,"in bytes");
	if(mm<m){
		printError(aD,"wrong value");
		return false;
	}

	xid->memorySizeMax = mm;

	xd->vecio = ippsMalloc_8u(mm);
	if(xd->vecio == NULL){
		printError(aD,"can not allocate memory for input/output files");
		return false;
	}
	return true;
}


bool form_files_names(allData *aD){
	char iform[30], oform[30];
	int sf, ss, sl, nm, ndmax, nv, nd;
	int ef, es, el;
	xInputData *xid;
	xOutputData *xod;

	sf = aD->row_first;
	sl = aD->row_last;
	ss = aD->row_step;
	nm = (sl-sf)/ss+1;
	nv = sl;
	ndmax = 1;
	while(nv >9){
		ndmax++;
		nv/=10;
	}

	xid = aD->hms->fbp->inputData;
	xod = aD->hms->fbp->outputData;

	nd = xid->NOD;
	
	if(nd>9){
		sprintf(aD->message,"number of digits for input files (%i, should be < 10)",nd);
		printError(aD);
		return false;
	}

	if(nd<ndmax){
		sprintf(aD->message,"number of digits for input files (%i), however the last image (%i) has %i digits",nd,sl,ndmax);
		printError(aD);
		return false;
	}

	if(strlen(xid->folder) < 1){
		sprintf(iform,"%%s%%%%0%ii%%s.%%s",nd);
		sprintf(aD->itemplate,iform,xid->prefix,xid->suffix,xid->extension);
	}else{
		sprintf(iform,"%%s/%%s%%%%0%ii%%s.%%s",nd);
		sprintf(aD->itemplate,iform,xid->folder,xid->prefix,xid->suffix,xid->extension);
	}
	
	ef = xod->fileFirst;
	es = xod->fileStep;
	el = ef+es*(nm-1);

	nv = el;
	ndmax = 1;
	while(nv >9){
		ndmax++;
		nv/=10;
	}

	nd = xod->NOD;
	
	if(nd>9){
		sprintf(aD->message,"number of digits for output files (%i, should be < 10)",nd);
		printError(aD);
		return false;
	}

	if(nd<ndmax){
		sprintf(aD->message,"number of digits for output files (%i), however the last image (%i) has %i digits",nd,el,ndmax);
		printError(aD);
		return false;
	}

	if(strlen(xod->folder) < 1){
		sprintf(oform,"%%s%%%%0%ii%%s.%%s",nd);
		sprintf(aD->otemplate,oform,xod->prefix,xod->suffix,xod->extension);
	}else{
		sprintf(oform,"%%s/%%s%%%%0%ii%%s.%%s",nd);
		sprintf(aD->otemplate,oform,xod->folder,xod->prefix,xod->suffix,xod->extension);
	}

	return true;
}





bool high_peaks_before(allData *aD){
	int nx, ny, nt, n;
	float jump;
	Ipp32f sv;
	xData *xd;
	xHighPeaks *hp;

	xd = aD->data;
	hp = aD->hms->fbp->preprocessing->highPeaksBefore;

	nx = aD->nx;
	ny = aD->ny;
	nt = nx * ny;
	
	printTagStart(aD,"HighIntensityPeaksBefore","High intensity peaks (before ring artefact suppression)");
	if(hp->type == HIGH_PEAKS_YES){
		printTag(aD,"Type","Yes","No, Yes");
	}else{
		printTag(aD,"Type","No","No, Yes");
	}
	if(hp->type == HIGH_PEAKS_NO){
		printTagEnd(aD);
		return true;
	}

	jump = hp->jump;

	printTag(aD,"Jump",jump,"jump between neighbouring pixels");
	if(jump<0){
		printError(aD,"should be a non-negative value");
		return false;
	}

	n = hp->numberPixels;
	printTag(aD,"Number",n,"number of neighboring pixels");
	if(n<1 || n> ny-1){
		printError(aD,"wrong number");
		return false;
	}
	if(aD->isFirstSlice){
		xd->vecp1 = ippsMalloc_32f(nt);
		xd->vecp2 = ippsMalloc_32f(nt);
		if(xd->vecp1 == NULL || xd->vecp2 == NULL){
			printError(aD,"can not allocate memory (high peaks)");
			return false;
		}
	}

	printTagStart(aD,"Steps");
	
	for(int k=1;k<=n;k++){
		ippsSub_32f(xd->veci+nx,xd->veci,xd->vecp1,(ny-1)*nx);
		ippsThreshold_LTValGTVal_32f(xd->vecp1, xd->vecp2, (ny-1)*nx, -jump, 1.0, -jump, 0.0);
		ippsThreshold_LT_32f_I(xd->vecp1,(ny-1)*nx,-jump);
		
		ippsSum_32f(xd->vecp2, (ny-1)*nx,&sv,ippAlgHintNone);
		printTag(aD,"Iteration",k);
		printTag(aD,"Direction","Top-to-bottom");
		printTag(aD,"NumberOfPixelsToBeCorrected",(int)sv);
		printTag(aD,"PercentOfPixelsToBeCorrected",float(100.0*sv/float((ny-1)*nx)));
		//fprintf(aD->fp,"Iteration #%i, Number of pixels to be corrected (top-to-bottom): %i (%f %%).\n",k,(int)sv,100.0*sv/float((ny-1)*nx));
		
		for(int i=0;i<ny-1;i++){
			ippsAdd_32f(xd->veci+(i+1)*nx,xd->vecp1+i*nx,xd->veci+i*nx,nx);
		}

		ippsSub_32f(xd->veci,xd->veci+nx,xd->vecp1,(ny-1)*nx);
		ippsThreshold_LTValGTVal_32f(xd->vecp1, xd->vecp2, (ny-1)*nx, -jump, 1.0, -jump, 0.0);
		ippsThreshold_LT_32f_I(xd->vecp1,(ny-1)*nx,-jump);
		
		ippsSum_32f(xd->vecp2, (ny-1)*nx,&sv,ippAlgHintNone);
		printTag(aD,"Direction","Bottom-to-top");
		printTag(aD,"NumberOfPixelsToBeCorrected",(int)sv);
		printTag(aD,"PercentOfPixelsToBeCorrected",float(100.0*sv/float((ny-1)*nx)));
		//fprintf(aD->fp,"Iteration #%i, Number of pixels to be corrected (top-to-bottom): %i (%f %%).\n",k,(int)sv,100.0*sv/float((ny-1)*nx));

		for(int i=ny-2;i>=0;i--){
			ippsAdd_32f(xd->veci+i*nx,xd->vecp1+i*nx,xd->veci+(i+1)*nx,nx);
		}
	}
	printTagEnd(aD);

	printTagEnd(aD);
	return true;
}


bool intensity_norm(allData *aD){
	int nx, ny, ntype, zl, zr, cl, cr;
	Ipp32f s, s0;
	xIntensity *xi;
	xData *xd;

	printTagStart(aD,"IntensityNormalisation");
	xd = aD->data;
	xi = aD->hms->fbp->preprocessing->intensity;
	nx = aD->nx;
	ny = aD->ny;
	ntype = xi->type;
	zl = xi->zeroLeft;
	zr = xi->zeroRight;
	cl = xi->columnLeft;
	cr = xi->columnRight;

	
	s0 = 0.0f;

	switch(ntype){
		case INTENSITY_NO:
			printTag(aD,"Type","No","no intensity normalisation");
			break;
		case INTENSITY_ROW:
			printTag(aD,"Type","Row","averaged optical path is the same for each row");
			zl = 0;
			zr = nx-1;
			ippsMean_32f(xd->veci+zl,zr-zl+1,&s0,ippAlgHintNone);
			break;
		case INTENSITY_COLUMN:
			printTag(aD,"Type","Column","averaged optical path is the same between the given columns");
			printTag(aD,"ColumnLeft",cl);
			printTag(aD,"ColumnRight",cr);
			if(cl<0 || cr >nx-1 || cl> cr){
				printError(aD,"wrong columns");
				return false;
			}
			zl = cl;
			zr = cr;
			ippsMean_32f(xd->veci+zl,zr-zl+1,&s0,ippAlgHintNone);
			break;
		case INTENSITY_ZERO:
			printTag(aD,"Type","Zero","optical path is zero between the given columns");
			printTag(aD,"ColumnLeft",zl);
			printTag(aD,"ColumnRight",zr);
			if(zl<0 || zr >nx-1 || zl> zr){
				printError(aD,"wrong columns");
				return false;
			}
			break;
		default:
		break;
	}

	if(ntype != INTENSITY_NO){
		for(int i=0;i<ny;i++){
			ippsMean_32f(xd->veci+i*nx+zl,zr-zl+1,&s,ippAlgHintNone);
			ippsSubC_32f_I(s-s0,xd->veci+i*nx,nx);
		}
	}

	printTagEnd(aD);
	return true;
}



bool find_missed_projections(allData *aD){
	int sl, mp;
	Ipp16u u;
	int v;
	bool is1, is2;
	xTransform *xt;
	xData *xd;

	xt = aD->hms->fbp->transform;
	xd = aD->data;

	sl = (int)strlen(xt->missedProjections);
	mp = 0;
	if(xt->missedProjectionsType != MISSED_PROJECTIONS_NO){
		if(sl >0){
			is2 = false;
			for(int i=0;i<sl;i++){
				v = (int)(*((Ipp8u*)(xt->missedProjections)+i)) - 48; 
				is1 = (v>=0 && v< 10);
				if(!is2 && is1)	mp++;
				is2 = is1;
			}
		}
	}

	aD->miss_proj = mp;
	xd->vmiss = ippsMalloc_16u(mp+1);
	if(xd->vmiss == NULL){
		printError(aD,"can not allocate memory (missed projections)");
		return false;
	}
	is1 = true;
	if(mp>0){	
		u = 0;
		is2 = false;
		mp = 0;
		for(int i=0;i<sl;i++){
			v = (int)(*((Ipp8u*)(xt->missedProjections)+i)) - 48; 
			is1 = (v>=0 && v< 10);
			if(!is2 && is1)	mp++;
			if(is1)	u = Ipp16u (10*u+ (Ipp16u)v);
			if(!is1 && is2){
				xd->vmiss[mp-1] = u;
				u = 0;
			}
			is2 = is1;
		}
		if(is1)	xd->vmiss[mp-1] = u;

		do{
			is1 = true;
			for(int i=1;i<mp;i++){
				if(xd->vmiss[i-1]<= xd->vmiss[i])continue;
				is1 = false;
				u = xd->vmiss[i-1];
				xd->vmiss[i-1] = xd->vmiss[i];
				xd->vmiss[i] = u;
			}
		}while(!is1);

		for(int i=mp-1;i>0;i--){
			if(xd->vmiss[i] == xd->vmiss[i-1]){
				for(int j= i;j<mp-1;j++){
					xd->vmiss[j] = xd->vmiss[j+1];
				}
				mp--;
			}
		}
		printTagStart(aD,"MissedProjections");

		for(int i=0;i<mp;i++){
			printTag(aD,"Projection",xd->vmiss[i]);
		}
		printTagEnd(aD);
		fprintf(aD->fp,"\n");
		aD->miss_proj = mp;

		if(xd->vmiss[mp-1]>= aD->ny){
			sprintf(aD->message,"wrong number of projection (%i), should be < %u",xd->vmiss[mp-1],aD->ny);
			printError(aD);
			return false;
		}
	}
	xd->vmiss[mp] = Ipp16u(aD->ny);
	return true;
}
bool check_crop(allData *aD){
	xTransform *xt;
	xt = aD->hms->fbp->transform;
	printTag(aD,"CropLeft",xt->cropLeft);
	if(xt->cropLeft>= aD->nx){
		printError(aD,"wrong number");
		return false;
	}
	printTag(aD,"CropRight",xt->cropRight);
	if(xt->cropRight>= aD->nx){
		printError(aD,"wrong number");
		return false;
	}

	if(xt->cropLeft + xt->cropRight >= aD->nx){
		sprintf(aD->message,"crop left + crop right (%u) >= image width (%u)",xt->cropLeft + xt->cropRight, aD->nx);
		printError(aD);
		return false;
	}


	printTag(aD,"CropTop",xt->cropTop);
	if(xt->cropTop>= aD->ny){
		printError(aD,"wrong number");
		return false;
	}
	printTag(aD,"CropBottom",xt->cropBottom);
	if(xt->cropBottom>= aD->ny){
		printError(aD,"wrong number");
		return false;
	}

	if(xt->cropTop + xt->cropBottom >= aD->ny){
		sprintf(aD->message,"crop top + crop bottom (%u) >= image height (%u)",xt->cropTop + xt->cropBottom, aD->ny);
		printError(aD);
		return false;
	}
	return true;
}


bool check_roi(allData *aD){
	int nx;
	xRoi *xr;
	xBackprojection *xb;
	xInputData *xi;
	xi = aD->hms->fbp->inputData;
	xb = aD->hms->fbp->backprojection;
	xr = xb->roi;
	nx = aD->nx;

	printTagStart2(aD,"ROI");
	if(xr->type == ROI_STANDARD){
		printTag(aD,"Type","Standard");
		xr->xmin = 0.0;
		xr->ymin = 0.0;
		if(xi->shape == SHAPE_POINT){
			xr->xmax = float(nx-1);
			xr->ymax = float(nx-1);
		}else{
			xr->xmax = float(nx);
			xr->ymax = float(nx);
		}
	}else if(xr->type == ROI_RECTANGLE){
		printTag(aD,"Type","Rectangle");

		if(xr->xmin < 0.0f){
			sprintf(aD->message,"x_min (%f), should be >= 0",xr->xmin);
			printWarning(aD);
			xr->xmin = 0.0f;
		}
		if(xr->ymin < 0.0f){
			sprintf(aD->message,"y_min (%f), should be >= 0",xr->ymin);
			printWarning(aD);
			xr->ymin = 0.0f;
		}
		if(xr->xmax > float(nx-1) && xi->shape == SHAPE_POINT){
			sprintf(aD->message,"x_max (%f), should be <= %i",xr->xmax,nx-1);
			printWarning(aD);
			xr->xmax = float(nx-1);
		}
		if(xr->xmax > float(nx) && xi->shape == SHAPE_PIXEL){
			sprintf(aD->message,"x_max (%f), should be <= %i",xr->xmax,nx);
			printWarning(aD);
			xr->xmax = float(nx);
		}
		if(xr->ymax > float(nx-1) && xi->shape == SHAPE_POINT){
			sprintf(aD->message,"ROI y_max (%f), should be <= %i",xr->ymax,nx-1);
			printWarning(aD);
			xr->ymax = float(nx-1);
		}
		if(xr->ymax > float(nx) && xi->shape == SHAPE_PIXEL){
			sprintf(aD->message,"y_max (%f), should be <= %i",xr->ymax,nx);
			printWarning(aD);
			xr->ymax = float(nx);
		}
		if(xr->xmax - xr->xmin < 1.0f){
			sprintf(aD->message,"x_min (%f), x_max (%f)",xr->xmin,xr->xmax);
			printError(aD);
			return false;
		}
		if(xr->ymax - xr->ymin < 1.0f){
			sprintf(aD->message,"y_min (%f), y_max (%f)",xr->ymin,xr->ymax);
			printError(aD);
			return false;
		}
	}
	printTag(aD,"Xmin",xr->xmin);
	printTag(aD,"Xmax",xr->xmax);
	printTag(aD,"Ymin",xr->ymin);
	printTag(aD,"Ymax",xr->ymax);

	printTagEnd2(aD);

	return true;
}



bool check_roi_size(allData *aD){
	int onx, ony;
//	xTransform *xt;
	xRoi *xr;
	gpu_info *gi;
	float w, h, ps;

	gi = aD->gi;

	ps = 1.0f;
	w = gi->xmax - gi->xmin;
	h = gi->ymax - gi->ymin;

//	xt = aD->hms->fbp->transform;
	xr = aD->hms->fbp->backprojection->roi;

	printTagStart2(aD,"ROIsize");

	printTag(aD,"WidthFloat",w);
	printTag(aD,"HeightFloat",h);

	gi->origWidth = w;
	gi->origHeight = h;

	if(xr->outputWidthType == OUTPUT_WIDTH_STANDARD){
		printTag(aD,"WidthType","Standard");
		ps = 1.0f;
		onx = int(w);
		ony = int(h);
	}else{
		printTag(aD,"WidthType","Given");
		ps = w/float(xr->outputWidth);
		onx = xr->outputWidth;
		ony = (int)(h*float(onx)/w);
	}
	gi->outputPixelSize = ps;
	printTag(aD,"OutputPixelSize",ps);
	printTag(aD,"WidthInt",onx);
	printTag(aD,"HeightInt",ony);

	printTagEnd2(aD);
	gi->mxi = onx;
	gi->myi = ony;
	return true;
}






bool printSettingsXml_FBP(allData *aD){
	printTagStart2(aD,"HMxml");
	printTagStart2(aD,"FBP");

	printTag(aD,"GPUDeviceNumber",aD->hms->fbp->GPUDeviceNumber);

	printTagStart2(aD,"BeamlineUser");
	printTagEnd2(aD);

	printTag(aD,"LogFile",aD->hms->fbp->logFile);

	xInputData *xi = aD->hms->fbp->inputData;
	printTagStart2(aD,"InputData");
	printTag(aD,"Folder",xi->folder);
	printTag(aD,"Prefix",xi->prefix);
	printTag(aD,"Suffix",xi->suffix);
	printTag(aD,"Extension",xi->extension);
	printTag(aD,"NOD",xi->NOD,"Number of digits");
	printTag(aD,"MemorySizeMax",xi->memorySizeFMax,"Maximal size of a file (MB)");
	printTag(aD,"MemorySizeMin",xi->memorySizeMin,"Minimasl size of a file (bytes)");

	switch(xi->orientation){
		case INPUT_ORIENTATION_SLICE:
			printTag(aD,"Orientation","Slice", "Slice, Projection");
			break;
		case INPUT_ORIENTATION_PROJECTION:
			printTag(aD,"Orientation","Projection", "Slice, Projection");
			break;
		default:
			break;
	}
		

	printTag(aD,"FileFirst",xi->fileFirst);
	printTag(aD,"FileLast",xi->fileLast);
	printTag(aD,"FileStep",xi->fileStep);
	printTag(aD,"ImageFirst",xi->imageFirst);
	printTag(aD,"ImageLast",xi->imageLast);
	printTag(aD,"ImageStep",xi->imageStep);

	printTagStart2(aD,"Raw");
	switch(xi->raw->type){
		case RAW_TYPE_YES:
			printTag(aD,"Type","Yes", "Yes, No");
			break;
		case RAW_TYPE_NO:
			printTag(aD,"Type","No", "Yes, No");
			break;
		default:
			break;
	}
	printTag(aD,"Bits",xi->raw->bits);
	switch(xi->raw->byteOrder){
		case BYTE_ORDER_BIG:
			printTag(aD,"ByteOrder","BE", "BE, LE");
			break;
		case BYTE_ORDER_LITTLE:
			printTag(aD,"ByteOrder","LE", "BE, LE");
			break;
		default:
			break;
	}
	printTag(aD,"Offset",xi->raw->offset);
	printTag(aD,"Gap",xi->raw->gap);
	printTag(aD,"Xlen",xi->raw->xlen);
	printTag(aD,"Ylen",xi->raw->ylen);
	printTag(aD,"Zlen",xi->raw->zlen);
	printTagEnd2(aD);
	
	printTag(aD,"FirstImageIndex",xi->firstImageIndex,"Index of the first image");
	printTag(aD,"ImagesPerFile",xi->imagesPerFile);

	switch(xi->restrictions){
		case INPUT_RESTRICTIONS_NO:
			printTag(aD,"Restrictions","No","No, Yes");
			break;
		case INPUT_RESTRICTIONS_YES:
			printTag(aD,"Restrictions","Yes","No, Yes");
			break;
		default:
			break;
	} 
	printTag(aD,"ValueMin",xi->valueMin);
	printTag(aD,"ValueMax",xi->valueMax);
	
	switch(xi->type){
		case INPUT_TYPE_INTENSITY:
			printTag(aD,"Type","Intensity","Intensity, Attenuation");
			break;
		case INPUT_TYPE_ATTENUATION:
			printTag(aD,"Type","Attenuation","Intensity, Attenuation");
			break;
		default:
			break;
	}

	switch(xi->shape){
		case SHAPE_POINT:
			printTag(aD,"Shape","Point","Point, Pixel");
			break;
		case SHAPE_PIXEL:
			printTag(aD,"Shape","Pixel","Point, Pixel");
			break;
		default:
			break;
	} 

	printTag(aD,"PixelParam",xi->pixelParam);
	
	printTagEnd2(aD);

	



	printTagStart2(aD,"FlatDarkFields");

	xFDField *xf = aD->hms->fbp->flatDarkFields->flatField;
	printTagStart2(aD,"FlatField");
	switch(xf->type){
		case FIELD_TYPE_USER:
			printTag(aD,"Type","User","User, Row");
			break;
		case FIELD_TYPE_ROW:
			printTag(aD,"Type","Row","User, Row");
			break;
		default:
			break;
	} 
	printTag(aD,"ValueBefore",xf->valueBefore);
	printTag(aD,"ValueAfter",xf->valueAfter);
	printTag(aD,"FileBefore",xf->fileBefore);
	printTag(aD,"FileAfter",xf->fileAfter);

	switch(xf->typeProfile){
		case FIELD_PROFILE_CONSTANT:
			printTag(aD,"ProfileType","Constant","Constant, Linear, Profile");
			break;
		case FIELD_PROFILE_LINEAR:
			printTag(aD,"ProfileType","Linear","Constant, Linear, Profile");
			break;
		case FIELD_PROFILE_PROFILE:
			printTag(aD,"ProfileType","Profile","Constant, Linear, Profile");
			break;
		default:
			break;
	} 

	printTag(aD,"FileProfile",xf->fileProfile);
	printTagEnd2(aD);

	xf = aD->hms->fbp->flatDarkFields->darkField;
	printTagStart2(aD,"DarkField");
	switch(xf->type){
		case FIELD_TYPE_USER:
			printTag(aD,"Type","User","User, Row");
			break;
		case FIELD_TYPE_ROW:
			printTag(aD,"Type","Row","User, Row");
			break;
		default:
			break;
	} 
	printTag(aD,"ValueBefore",xf->valueBefore);
	printTag(aD,"ValueAfter",xf->valueAfter);
	printTag(aD,"FileBefore",xf->fileBefore);
	printTag(aD,"FileAfter",xf->fileAfter);

	switch(xf->typeProfile){
		case FIELD_PROFILE_CONSTANT:
			printTag(aD,"ProfileType","Constant","Constant, Linear, Profile");
			break;
		case FIELD_PROFILE_LINEAR:
			printTag(aD,"ProfileType","Linear","Constant, Linear, Profile");
			break;
		case FIELD_PROFILE_PROFILE:
			printTag(aD,"ProfileType","Profile","Constant, Linear, Profile");
			break;
		default:
			break;
	} 

	printTag(aD,"FileProfile",xf->fileProfile);
	printTagEnd2(aD);
	printTagEnd2(aD);



	printTagStart2(aD,"Preprocessing");
	printTagEnd2(aD);

	printTagStart2(aD,"Transform");
	printTagEnd2(aD);

	printTagStart2(aD,"Backprojection");
	printTagEnd2(aD);

	
	xOutputData *xo = aD->hms->fbp->outputData;
	printTagStart2(aD,"OutputData");
	switch(xo->type){
		case OUTPUT_TYPE_SOLUTION:
			printTag(aD,"Type","Solution","Solution, FFCorrection, Preprocessed, Transformed");
			break;
		case OUTPUT_TYPE_FFCORRECTION:
			printTag(aD,"Type","FFCorrection","Solution, FFCorrection, Preprocessed, Transformed");
			break;
		case OUTPUT_TYPE_PREPROCESSED:
			printTag(aD,"Type","Preprocessed","Solution, FFCorrection, Preprocessed, Transformed");
			break;
		case OUTPUT_TYPE_TRANSFORMED:
			printTag(aD,"Type","Transformed","Solution, FFCorrection, Preprocessed, Transformed");
			break;
		default:
			break;
	} 

	switch(xo->state){
		case OUTPUT_STATE_ATTENUATION:
			printTag(aD,"State","Attenuation","Attenuation, Intensity");
			break;
		case OUTPUT_STATE_INTENSITY:
			printTag(aD,"State","Intensity","Attenuation, Intensity");
			break;
		default:
			break;
	} 

	printTag(aD,"Folder",xo->folder);
	printTag(aD,"Prefix",xo->prefix);
	printTag(aD,"Suffix",xo->suffix);
	printTag(aD,"Extension",xo->extension);
	printTag(aD,"NOD",xo->NOD);
	printTag(aD,"FileFirst",xo->fileFirst);
	printTag(aD,"FileStep",xo->fileStep);
	switch(xo->bitsType){
		case OUTPUT_BITS_INPUT:
			printTag(aD,"BitsType","Input","Input, Given");
			break;
		case OUTPUT_BITS_GIVEN:
			printTag(aD,"BitsType","Given","Input, Given");
			break;
		default:
			break;
	} 
	printTag(aD,"Bits",xo->bits);
	switch(xo->restrictions){
		case OUTPUT_RESTRICTIONS_NO:
			printTag(aD,"Restrictions","No","No, Yes");
			break;
		case OUTPUT_RESTRICTIONS_YES:
			printTag(aD,"Restrictions","Yes","No, Yes");
			break;
		default:
			break;
	} 
	printTag(aD,"ValueMin",xo->valueMin);
	printTag(aD,"ValueMax",xo->valueMax);
	switch(xo->shape){
		case SHAPE_POINT:
			printTag(aD,"Shape","Point","Point, Pixel, Input");
			break;
		case SHAPE_PIXEL:
			printTag(aD,"Shape","Pixel","Point, Pixel, Input");
			break;
		case SHAPE_INPUT:
			printTag(aD,"Shape","Input","Point, Pixel, Input");
			break;
		default:
			break;
	} 
	printTagEnd2(aD);

	printTagEnd2(aD);
	printTagEnd2(aD);

	return true;
}






bool find_params(allData *aD){
	int xbl, ybl;
	float xd, yd, diam;
	float xc;
	float roi_orig_x, roi_orig_y;
	float angle_roi;
	float dm;
	int idm, uu;
	int ux, uy;
	int suy, sux, op;
	float sdm;
	float wr, hr;

	gpu_info *gi;
	xBackprojection *xb;
	xRoi *xr;
	
	xr = aD->hms->fbp->backprojection->roi;
	xb = aD->hms->fbp->backprojection;


	gi = aD->gi;
	int_ratio(gi->mxi,gi->wo,&(gi->mxo));
	int_ratio(gi->myi,gi->ho,&(gi->myo));
	printTagStart2(aD,"GPUThreads");
	printTag(aD,"XThreadsPerBlock",gi->wo);
	printTag(aD,"YThreadsPerBlock",gi->ho);
	printTag(aD,"WidthBefore",gi->mxi);
	printTag(aD,"Width",gi->mxo);
	printTag(aD,"HeightBefore",gi->myi);
	printTag(aD,"Height",gi->myo);
	xbl = gi->mxo/gi->wo;
	ybl = gi->myo/gi->ho;
	printTag(aD,"XBlocks",xbl);
	printTag(aD,"YBlocks",ybl);
	printTag(aD,"TotalBlocks",xbl*ybl);

	xd = gi->outputPixelSize * float(gi->wo);
	yd = gi->outputPixelSize * float(gi->ho);
	diam =  sqrt(xd*xd + yd*yd);
	gi->outputBlockDiameter = diam;
	printTag(aD,"OutputBlockDiameter",diam);

	wr = gi->mxo * gi->outputPixelSize;
	hr = gi->myo * gi->outputPixelSize;

	gi->outputROIRadius = 0.5f*sqrt(wr*wr + hr*hr);
	printTag(aD,"OutputROIRadius",gi->outputROIRadius);

	
	
	uu = gi->wo*gi->ho;
	ux = uu;
	uy = 1;

	idm = (int)ceil(diam)+2;

	if(idm > uu){
		sprintf(aD->message,"the diameter of the block is too large (%f)",diam);
		printError(aD);
		return false;
	}
	
	op = 0;
	dm = 2;
	sux = 2;
	suy = 1;
	printTagStart2(aD,"InputChunkOptions");
	for(int i=1; i<=uu; i++){
		if(uu%i != 0)continue;
		sprintf(aD->message,"Option%i",op+1);
		printTagStart2(aD,aD->message);
		ux = uu/i;
		printTag(aD,"Width",ux);
		uy = i;
		printTag(aD,"Height",uy);
		dm = diam;// + gi->outputROIRadius*sinf(float(uy-1)*gi->rotAngleStep);
		printTag(aD,"Diameter",dm);
		idm = (int)ceil(dm);
		printTag(aD,"DiameterInteger",idm);
		if(idm + 3 < ux){
			printTag(aD,"Result","true");
			sdm = dm;
			sux = ux;
			suy = uy;
		}else{
			printTag(aD,"Result","false");
			printTagEnd2(aD);//Option
			break;
		}
		op++;
		printTagEnd2(aD);//Option
	}
	printTagEnd2(aD);//InputChunkOptions

	ux = sux;
	uy = suy;
	sdm = dm;
	
	
	gi->blockRadius = 0.5f*sdm;

	printTag(aD,"BlockRadius",gi->blockRadius);
	printTag(aD,"InputBlockWidth",ux);
	printTag(aD,"InputBlockHeight",uy);

	gi->ux = ux;
	gi->uy = uy;

	//printf("New ux: %i, uy: %i\n",ux,uy);
	//printf("New ux: %i, uy: %i\n",aD->gi->ux,aD->gi->uy);


	xc = xb->imageCentre;
	printTag(aD,"ImageCentre",xc);

	roi_orig_x = gi->xmin + gi->outputPixelSize *0.5f*float(gi->mxo);
	roi_orig_y = gi->ymin + gi->outputPixelSize *0.5f*float(gi->myo);

	printTag(aD,"ROIOriginalXNoRotation",roi_orig_x);
	printTag(aD,"ROIOriginalYNoRotation",roi_orig_y);

	float x1, y1;
	x1 = roi_orig_x - xc;
	y1 = roi_orig_y - xc;

	angle_roi = xr->angle * pif/180.0f;

	roi_orig_x = xc + (x1*cosf(angle_roi)) - (y1*sinf(angle_roi));
	roi_orig_y = xc + (x1*sinf(angle_roi)) + (y1*cosf(angle_roi));

	printTag(aD,"ROIOriginalX",roi_orig_x);
	printTag(aD,"ROIOriginalY",roi_orig_y);

	//ddx = 0.5f*(gi->outputPixelSize * float(gi->mxo) - gi->origWidth);
	//ddy = 0.5f*(gi->outputPixelSize * float(gi->myo) - gi->origHeight);
	

//	roi_new_x = roi_orig_x + ddx*cos(angle_roi) - ddy*sin(angle_roi);
//	roi_new_y = roi_orig_y + ddx*sin(angle_roi) + ddy*cos(angle_roi);

//	gi->newX = roi_new_x;
//	gi->newY = roi_new_y;
	gi->angle = angle_roi;
	gi->newWidth = gi->outputPixelSize * float(gi->mxo);
	gi->newHeight = gi->outputPixelSize * float(gi->myo);

//	printTag(aD,"ROINewX",roi_new_x);
//	printTag(aD,"ROINewY",roi_new_y);

	float roiR, roiA, roiRad;
	x1 = roi_orig_x - xc;
	y1 = roi_orig_y - xc;
	roiRad = 0.5f*sqrt(gi->newWidth * gi->newWidth + gi->newHeight * gi->newHeight);

	printTag(aD,"ROINewXrel",x1);
	printTag(aD,"ROINewYrel",y1);
	roiR = sqrt(x1*x1+y1*y1);
	roiA = atan2(y1,x1);
	gi->roiR = roiR;
	gi->roiA = roiA;
	printTag(aD,"ROINewRrel",roiR);
	printTag(aD,"ROINewArel",180.0f*roiA/pif);
	printTag(aD,"ROIRad",roiRad);

	int x_min1, x_max1, x_min2, x_max2;
	x_min1 = (int)floor(xc - roiRad)-1;
	x_max1 = (int)ceil(xc + roiRad)+1;
	printTag(aD,"ROIXMin",x_min1);
	printTag(aD,"ROIXMax",x_max1);

	xTransform *xt;
	xt = aD->hms->fbp->transform;

	x_min2 = xt->cropLeft;
	x_max2 = aD->nx - xt->cropRight;
	if(aD->hms->fbp->inputData->shape == SHAPE_PIXEL) x_max2++;
	printTag(aD,"OrigXMin",x_min2);
	printTag(aD,"OrigXMax",x_max2);

	gi->x_min = __min(x_min1,x_min2);
	gi->x_max = __max(x_max1,x_max2);

	printTag(aD,"FinalXMin",gi->x_min);
	printTag(aD,"FinalXMax",gi->x_max);




	int extpower, next, origw;

	origw = gi->x_max - gi->x_min;
	gi->origw = origw;
	printTag(aD,"RequiredWidth",origw);
	find_power(origw,&next,&extpower);
	
	gi->fpower = extpower;
	gi->flength = next;

	printTag(aD,"WidthAfterExtrapolation",next);
	printTag(aD,"PowerOf2WidthAfterExtrapolation",extpower);
	xt->extrapolationWidth = next;




	int v, vH;
	v = gi->wo * gi->ho;
	nom(v,gi->warpSize, &vH);
	vH/=2;
	printTag(aD,"SinCosChunk",vH);
	gi->vertH = vH;

	gi->nuc = vH;
	if(uy>0){
		gi->nuc = vH/uy;
	}
	printTag(aD,"BlocksInSinCosChunk",gi->nuc);

	printTagEnd2(aD);

	return true;
}

