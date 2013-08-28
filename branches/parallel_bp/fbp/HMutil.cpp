/* ----------------------------------------------------------------------
 * HMutil.cpp()
 * ----------------------------------------------------------------------
 * 
 * ----------------------------------------------------------------------
 */
#include "defs.h"


bool readInputVolume(MatPar *matp, PolarToCart *PC, allData *aD, OtherParam *other){
	int nx, ny, nz, nt, ew, ny_old, r;
	nx = aD->nx;
	nz = aD->nz;
	ny = matp->y_maxi - matp->y_mini+1;
	nt = nx*ny*nz;
	ew = other->eWidth;
	ny_old = aD->ny;
	r = matp->y_mini - aD->hms->af->rowFirst;

	printf("ew: %i\n", ew);
	
	matp->ivol = ippsMalloc_32f(nt);
	matp->rec = ippsMalloc_32f(PC->OutputWidth  * PC->OutputHeight);
	matp->vsect = ippsMalloc_32fc(ew * ny);
	matp->fsect = ippsMalloc_32fc(ew * ny);

	if(matp->ivol == NULL || matp->rec == NULL || matp->vsect == NULL || matp->fsect == NULL){
		fprintf(aD->fp, "Error: can not allocate memory (volume).\n");
		return false;
	}

	if(aD->hms->af->clockwiseRotation == CLOCKWISE_ROTATION_YES){
		for(int i=0;i<nz;i++){
			ippsCopy_32f(aD->data->veci+nx*ny_old*i+ r*nx,matp->ivol+nx*ny*i,nx*ny);
		}
	}else{
		for(int i=0;i<nz;i++){
			ippsCopy_32f(aD->data->veci+nx*ny_old*i+ r*nx,matp->ivol+nx*ny*(nz-1-i),nx*ny);
		}
	}
	
	return true;
}

bool fillMatPar(MatPar *matp, PolarToCart *PC, allData *aD, OtherParam *other, float row){
	xAllFDK *af = aD->hms->af;
	matp->na = other->num_images+1;
	matp->nr = PC->nrn;
	matp->rmin = PC->rmin;
	matp->rmax = PC->rmax;
	matp->ImageWidth = aD->nx;
	matp->ImageHeight = (int)af->originalImageHeight;
	matp->SourceToDetector = af->sourceToObject + af->detectorToObject;
	matp->PixelSize = af->pixelSize;
	matp->rowf = row;
	matp->ConeCentreShiftX = af->sourceXShift;
	matp->ConeCentreShiftY = af->sourceYShift;
	matp->AxisShift = af->axisXShift;
	matp->RestrMin = af->rowFirst;
	matp->RestrMax = af->rowFirst + aD->ny-1;
	return true;
}


bool formRowFilter(allData *aD, OtherParam *other){
	other->vecFilterRowR = ippsMalloc_32f(other->eWidth);
	//Ipp32f *zero = ippsMalloc_32f(other->eWidth);
	//Ipp32f *vec = ippsMalloc_32f(other->eWidth);
	float omega = aD->hms->af->filterParam;
	Ipp32f ds =1.0f;
	formFilterRL2(other->vecFilterRowR,omega, ds, other->eWidth);
	//ippsRealToCplx_32f(vec, zero, other->vecFilterRowC, other->eWidth);
	//ippsFree(vec); vec = NULL;
	//ippsFree(zero); zero = NULL;
	return true;
}

bool findSlices(int *SliceFirst, int *SliceLast, int *SliceStep, int RestrFirst, int RestrLast){
	int sf, sl, ss, rf, rl;
	sf = *SliceFirst;
	sl = *SliceLast;
	ss = *SliceStep;
	rf = RestrFirst;
	rl = RestrLast;
	rf = __max(rf,sf);
	rl = __min(rl,sl);
	if(rf > rl)return false;
	if(ss <1) return false;
	while(sf < rf){
		sf+=ss;
	};
	sl = rl;
	if(sl<sf)return false;
	sl = sf;
	do{
		sl+=ss;
	}while(sl<=rl);
	sl-=ss;
	*SliceFirst = sf;
	*SliceLast = sl;
	return true;
}

bool fillPC(PolarToCart *PC, allData *aD, OtherParam *other){
	float xmin, xmax, ymin, ymax, nx;
	xAllFDK *af = aD->hms->af;
	xmin = af->ROIXmin;
	xmax = af->ROIXmax;
	ymin = af->ROIYmin;
	ymax = af->ROIYmax;
	nx = float(aD->nx);
	fprintf(aD->fp,"ROI: x_min (%f), x_max (%f), y_min (%f), y_max (%f).\n",xmin,xmax,ymin,ymax);

	if(xmax < xmin+10){
		fprintf(aD->fp,"Error: wrong numbers (x).\n");
		return false;
	}

	if(ymax < ymin+10){
		fprintf(aD->fp,"Error: wrong numbers (y).\n");
		return false;
	}

	if(xmin <0 || xmax > float(nx-1)){
		fprintf(aD->fp,"Error: wrong numbers (x, edges).\n");
		return false;
	}

	if(ymin <0 || ymax > float(nx-1)){
		fprintf(aD->fp,"Error: wrong numbers (y, edges).\n");
		return false;
	}

	fprintf(aD->fp,"Rotation angle: %f (degrees).\n",af->rotationAngle);

	float ow;
	int outw;

	if(af->outputWidthType == OUTPUT_WIDTH_STANDARD){
		ow = xmax - xmin;
		outw = int(ow);
	}else{
		outw = af->outputWidth;
	}

	fprintf(aD->fp,"Output width: %i.\n",outw);

	fprintf(aD->fp,"Number of circles: %i.\n",af->numberOfCircles);
	if(af->numberOfCircles < 10){
		fprintf(aD->fp,"Error: wrong value.\n");
		return false;
	}
	
	PC->xmin = xmin;
	PC->xmax = xmax;
	PC->ymin = ymin;
	PC->ymax = ymax;
	PC->RotAngle = af->rotationAngle;
	PC->OutputWidth = outw;
	PC->num_images = other->num_images;
	PC->NumberOfCircles = af->numberOfCircles;
	PC->ImageWidth = aD->nx;
	PC->l_a = other->l_a;
	PC->l_r = other->l_r;
	PC->w_a = other->w_a;
	PC->w_r = other->w_r;
	return true;
}

bool ifFolderExists(allData *aD, char *fname){
	STAT_TYPE buf;
	int result;

	sprintf(aD->message,"Checking if folder \"%s\" exists",fname);
	printInfo(aD);
	
	result = STAT_FUNC(fname, &buf );
  
	if( result != 0 )
	{
		printWarning(aD, "Problem getting information about the folder");
		switch (errno)
		{
		case ENOENT:
			printWarning(aD,"Folder does not exist");
			break;
		case EINVAL:
			printWarning(aD,"Invalid parameter to _stat");
			break;
		default:
			printWarning(aD,"Unexpected error in _stat");
			break;
		}
		return false;
	}

	printInfo(aD,"The folder exists");
	return true;
}

bool createFolder(allData *aD, char *fname){
	int result;
	sprintf(aD->message,"Trying to create folder \"%s\"",fname);
	printInfo(aD);
	
	result = MKDIR_FUNC(fname);
  
	if( result != 0 )
	{
		printWarning(aD,"Can not create the folder");
		switch (errno)
		{
		case EEXIST:
			printWarning(aD,"Its name is the name of an existing file, directory, or device");
			break;
		case ENOENT:
			printWarning(aD,"The path has not been found");
			break;
		default:
			printWarning(aD,"Unexpected error");
			break;
		}
		return false;
	}

	printInfo(aD, "The folder has been created");
	return true;
}

bool createFolderCycle(allData *aD, char *fname){
	char tname[MAX_FOLDER];
	int i, m;

	i = 0;
	m = 2;
	
	if(fname[0] == '/' && fname[1] == '/') m =0;
	for(;;){
		if(m > 0 && (fname[i] == '/' || fname[i] == '\0')){
			if(m != 1){
				tname[i] = '\0';
			}else{
				tname[i] = '/';
				tname[i+1] = '\0';
			}
				
			if(!ifFolderExists(aD,tname)){
				if(!createFolder(aD,tname)){
					sprintf(aD->message, "can not create the folder \"%s\"",tname);
					printWarning(aD);
					return false;
				};
			}
			m++;
		}
		if(m == 0 && i>1 && fname[i] == '/') m = 1;
		if(fname[i] == ':' && (fname[i+1] == '/')){
			tname[i] = fname[i];
			tname[i+1] = fname[i+1];
			i++;
			tname[i+1] = '\0';
			if(!ifFolderExists(aD,tname)){
				if(!createFolder(aD,tname)){
					sprintf(aD->message, "can not create the folder \"%s\"",tname);
					printWarning(aD);
					return false;
				};
			}
		}
		if(fname[i] == '\0')return true;
		
		tname[i] = fname[i];
		i++;
	}
}


