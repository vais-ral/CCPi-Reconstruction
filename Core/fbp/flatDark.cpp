/* ----------------------------------------------------------------------
 * flatDark.cpp()
 * ----------------------------------------------------------------------
 * Reading flat/dark fields
 * Flat field correction
 * ----------------------------------------------------------------------
 */

#include "defs.h"


bool ff_correction(allData *aD){
	printTagStart(aD,"FlatFieldCorrection");
	//fprintf(aD->fp,"\n>> Flat field correction.\n");
	//printf("\n>> Flat field correction.\n");

	xFDField *xfdf;
	float w;
	xData *xd;
	xd = aD->data;
	
	
	xfdf = aD->hms->fbp->flatDarkFields->darkField;
	
	ippsCopy_32f(xd->matd1+aD->row*aD->nx,xd->vect1,aD->nx);
	ippsCopy_32f(xd->matd2+aD->row*aD->nx,xd->vect2,aD->nx);

	switch(xfdf->typeProfile){
		case FIELD_PROFILE_CONSTANT:
			for(unsigned int i=0;i<aD->ny;i++){
				ippsCopy_32f(xd->vect1,xd->vecd+i*aD->nx,aD->nx);
			}
			break;
		case FIELD_PROFILE_LINEAR:
			for(unsigned int i=0;i<aD->ny;i++){
				if(aD->ny <2){
					w = 0.0;
				}else{
					w = float(i)/float(aD->ny-1);
				}
				ippsMulC_32f(xd->vect2,w,xd->vecd+i*aD->nx,aD->nx);
				ippsAddProductC_32f(xd->vect2,1.0f-w,xd->vecd+i*aD->nx,aD->nx);
			}
			break;
		case FIELD_PROFILE_PROFILE:
			for(unsigned int i=0;i<aD->ny;i++){
				ippsMulC_32f(xd->vect1,xd->profd[aD->row*aD->ny+i],xd->vecd+i*aD->nx,aD->nx);
			}
			break;
		default:
			fprintf(aD->fp,"Something is wrong in ff correction.\n");
			break;
	}
	

	xfdf = aD->hms->fbp->flatDarkFields->flatField;
	
	ippsCopy_32f(xd->matf1+aD->row*aD->nx,xd->vect1,aD->nx);
	ippsCopy_32f(xd->matf2+aD->row*aD->nx,xd->vect2,aD->nx);

	switch(xfdf->typeProfile){
		case FIELD_PROFILE_CONSTANT:
			for(unsigned int i=0;i<aD->ny;i++){
				ippsCopy_32f(xd->vect1,xd->vecf+i*aD->nx,aD->nx);
			}
			break;
		case FIELD_PROFILE_LINEAR:
			for(unsigned int i=0;i<aD->ny;i++){
				if(aD->ny <2){
					w = 0.0;
				}else{
					w = float(i)/float(aD->ny-1);
				}
				ippsMulC_32f(xd->vect2,w,xd->vecf+i*aD->nx,aD->nx);
				ippsAddProductC_32f(xd->vect2,1.0f-w,xd->vecf+i*aD->nx,aD->nx);
			}
			break;
		case FIELD_PROFILE_PROFILE:
			for(unsigned int i=0;i<aD->ny;i++){
				ippsMulC_32f(xd->vect1,xd->proff[aD->row*aD->ny+i],xd->vecf+i*aD->nx,aD->nx);
			}
			break;
		default:
			fprintf(aD->fp,"Something is wrong in ff correction.\n");
			break;
	}

	

	ippsSub_32f_I(xd->vecd,xd->vecf,aD->nt);
	ippsSub_32f_I(xd->vecd,xd->veci,aD->nt);
	ippsThreshold_LT_32f_I(xd->vecf,aD->nt,0.1f);
	ippsThreshold_LT_32f_I(xd->veci,aD->nt,0.1f);
	ippsDiv_32f_I(xd->vecf,xd->veci,aD->nt);
	ippsLn_32f_I(xd->veci,aD->nt);
	ippsMulC_32f_I(-1.0,xd->veci,aD->nt);
	

	printTagEnd(aD);
	//printf("<< Flat field correction.\n\n");

	return true;

}





bool read_field(allData *aD, int df){
	char sname[20], uname[20];
	xFDField *xf;
	unsigned int nx1, ny1, nb1;
	unsigned int nx2, ny2, nb2, nt;
	gtiff *ti1, *ti2, *ti3;
	Ipp32f *mat1, *mat2, *prof;
	unsigned int nxp, nyp, nbp, ntp;
	xData *xd;

	xd = aD->data;

	if(df == READ_DARK){
		sprintf(sname,"dark");
		xf = aD->hms->fbp->flatDarkFields->darkField;
		printTagStart(aD,"DarkField");
	}else{
		sprintf(sname,"flat");
		xf = aD->hms->fbp->flatDarkFields->flatField;
		printTagStart(aD,"FlatField");
	}
	sprintf(uname,"%s",sname);
	toUpper(uname);
		
	switch(xf->type){
		case FIELD_TYPE_USER:
			printTag(aD,"Type","constant defined by the user");
			printInfo(aD,"the matrix will be set later");
			printTagEnd(aD);
			return true;
		case FIELD_TYPE_ROW:
			printTag(aD,"Type","row");
			break;
		default:
			sprintf(aD->message,"unknown type (%i)",xf->type);
			printError(aD);
			return false;
	}
	
	ti1 = new gtiff();
	ti2 = new gtiff();
	ti3 = new gtiff();
	if(!ti1->getInfo(xf->fileBefore,&nx1,&ny1,&nb1)){
		sprintf(aD->message,"can not read data from the %s field file \"%s\"",sname,xf->fileBefore);
		printError(aD);
		return false;
	}
		
	nt = nx1*ny1;

	if(df == READ_DARK){
		aD->dark_nx = nx1;
		aD->dark_ny = ny1;
	}else{
		aD->flat_nx = nx1;
		aD->flat_ny = ny1;
	}
	
	if(xf->typeProfile != FIELD_PROFILE_CONSTANT){
		if(!ti2->getInfo(xf->fileAfter,&nx2,&ny2,&nb2)){
			sprintf(aD->message,"can not read data from the %s field file \"%s\"",sname,xf->fileAfter);
			printError(aD);
			return false;
		}
				
		if(nx1 != nx2 || ny1 != ny2){
			sprintf(aD->message,"different width/height for the %s files",sname);
			printError(aD);
			return false;
		}
	}
	
	if(nt <= 0){
		printError(aD,"wrong size");
		return false;
	}

	if(df == READ_DARK){
		xd->matd1 = ippsMalloc_32f(nt);
		mat1 = xd->matd1;
		xd->matd2 = ippsMalloc_32f(nt);
		mat2 = xd->matd2;
	}else{
		xd->matf1 = ippsMalloc_32f(nt);
		mat1 = xd->matf1;
		xd->matf2 = ippsMalloc_32f(nt);
		mat2 = xd->matf2;
	}
	if(mat1 == NULL || mat2 == NULL){
		sprintf(aD->message,"can not allocate memory for the data from the %s field file(s)",sname);
		printError(aD);
		return false;
	}
	ti1->read_float(xf->fileBefore,mat1);
	delete ti1;

	if(xf->typeProfile != FIELD_PROFILE_CONSTANT){
		ti2->read_float(xf->fileAfter,mat2);
	}else{
		ippsCopy_32f(mat1,mat2,nt);
	}
	delete ti2;
	

	if(xf->typeProfile == FIELD_PROFILE_PROFILE){
		if(!ti3->getInfo(xf->fileProfile,&nxp,&nyp,&nbp)){
			sprintf(aD->message,"can not read data from the %s profile file \"%s\"",sname,xf->fileProfile);
			printError(aD);
			return false;
		}
		if(nbp != 32){
			sprintf(aD->message,"wrong number of bits (%u) in %s profile file \"%s\"",nbp,sname,xf->fileProfile);
			printError(aD);
			return false;
		}
		ntp = nxp * nyp;
		if(df == READ_DARK){
			xd->profd = ippsMalloc_32f(ntp);
			prof = xd->profd;
			aD->profd_nx = nxp;
			aD->profd_ny = nyp;
		}else{
			xd->proff = ippsMalloc_32f(ntp);
			prof = xd->proff;
			aD->proff_nx = nxp;
			aD->proff_ny = nyp;
		}
		if(prof == NULL){
			sprintf(aD->message,"can not allocate memory to read the %s profile file \"%s\"",sname,xf->fileProfile);
			printError(aD);
			return false;
		}
		ti3->read_float(xf->fileProfile,prof);
	}

	delete ti3;
	printTagEnd(aD);
	return true;
}


bool read_fdf_data(allData *aD){
	int nt;
	xData *xd;
	xFlatDarkFields * xfdf;

	xd = aD->data;

	printTagStart(aD,"ReadDarkFlat","Reading dark noise/flat field files");
			
	if (!read_field(aD, READ_DARK))return false;
	if (!read_field(aD, READ_FLAT))return false;
	
	xfdf = aD->hms->fbp->flatDarkFields;
	
	if(xfdf->darkField->type == FIELD_TYPE_ROW && xfdf->flatField->type == FIELD_TYPE_ROW){
		if(aD->dark_nx != aD->flat_nx){
			sprintf(aD->message,"different widths of dark (%u) and flat (%u) fields",aD->dark_nx,aD->flat_nx);
			printError(aD);
			return false;
		}
		if(aD->dark_ny != aD->flat_ny){
			sprintf(aD->message,"different heights of dark (%u) and flat (%u) fields",aD->dark_ny,aD->flat_ny);
			printError(aD);
			return false;
		}
	}else if(xfdf->darkField->type == FIELD_TYPE_ROW || xfdf->flatField->type == FIELD_TYPE_ROW){
		if(xfdf->darkField->type == FIELD_TYPE_ROW){
			nt = aD->dark_nx * aD->dark_ny;
			aD->flat_nx = aD->dark_nx;
			aD->flat_ny = aD->dark_ny;
			xd->matf1 = ippsMalloc_32f(nt);
			xd->matf2 = ippsMalloc_32f(nt);
			if(xd->matf1 == NULL || xd->matf2 == NULL){
				printError(aD,"can not allocate memory (flat field)");
				return false;
			}
			ippsSet_32f(xfdf->flatField->valueBefore,xd->matf1,nt);
			if(xfdf->flatField->typeProfile == FIELD_PROFILE_CONSTANT){
				ippsSet_32f(xfdf->flatField->valueBefore,xd->matf2,nt);
			}else{
				ippsSet_32f(xfdf->flatField->valueAfter,xd->matf2,nt);
			}
		}else{
			nt = aD->flat_nx * aD->flat_ny;
			aD->dark_nx = aD->flat_nx;
			aD->dark_ny = aD->flat_ny;
			xd->matd1 = ippsMalloc_32f(nt);
			xd->matd2 = ippsMalloc_32f(nt);
			if(xd->matd1 == NULL || xd->matd2 == NULL){
				printError(aD,"can not allocate memory (dark field)");
				return false;
			}
			ippsSet_32f(xfdf->darkField->valueBefore,xd->matd1,nt);
			
			if(xfdf->darkField->typeProfile == FIELD_PROFILE_CONSTANT){
				ippsSet_32f(xfdf->darkField->valueBefore,xd->matd2,nt);
			}else{
				ippsSet_32f(xfdf->darkField->valueAfter,xd->matd2,nt);
			}
		}
	}

	if(xfdf->darkField->type == FIELD_TYPE_ROW && xfdf->flatField->type == FIELD_TYPE_ROW){
		if(aD->dark_nx != aD->flat_nx){
			sprintf(aD->message,"different widths of dark (%u) and flat (%u) fields",aD->dark_nx,aD->flat_nx);
			printError(aD);
			return false;
		}
		if(aD->dark_ny != aD->flat_ny){
			sprintf(aD->message,"different heights of dark (%u) and flat (%u) fields",aD->dark_ny,aD->flat_ny);
			printError(aD);
			return false;
		}
	}
	
	printTagEnd(aD);
		
	return true;
}



