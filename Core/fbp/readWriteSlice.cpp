/* ----------------------------------------------------------------------
 * readWriteSlice.cpp()
 * ----------------------------------------------------------------------
 * Reading input slice (sinogram), writing output image
 * ----------------------------------------------------------------------
 */

#include "defs.h"
#include "time_stamp.h"

bool read_one_slice_fs(allData *aD){
	int nx, ny, nz, nb, nxz;
	float vmin, vmax;
	xFDField *xff, *xdf;
	xData *xd;
	xInputData *xid;

	xd = aD->data;
	xid = aD->hms->fbp->inputData;
	xdf = aD->hms->fbp->flatDarkFields->darkField;
	xff = aD->hms->fbp->flatDarkFields->flatField;
	
	nx = aD->temp_nx;
	ny = aD->temp_ny;
	nb = aD->temp_nb;
	
	
	if(nx <= 0){
		printError(aD,"wrong width");
		return false;
	}
	if(ny <= 0){
		printError(aD,"wrong height");
		return false;
	}
	
	if(!(nb == 8 || nb == 16 || nb == 32)){
		printError(aD,"wrong number of bits");
		return false;
	}

	aD->nx = nx;
	aD->ny = ny;
	aD->nb = nb;
	aD->nt = nx * ny;

	nz = aD->row_last + 1;
	nxz = nx * nz;
		
	if(xdf->type == FIELD_TYPE_USER && xff->type == FIELD_TYPE_USER){
		xd->matf1 = ippsMalloc_32f(nxz);
		xd->matf2 = ippsMalloc_32f(nxz);
		xd->matd1 = ippsMalloc_32f(nxz);
		xd->matd2 = ippsMalloc_32f(nxz);

		if(xd->matd1 == NULL || xd->matd2 == NULL || xd->matf1 == NULL || xd->matf2 == NULL){
			printError(aD,"can not allocate memory for flat/dark fields");
			return false;
		}
		aD->dark_nx = nx;
		aD->dark_ny = nz;
		aD->flat_nx = nx;
		aD->flat_ny = nz;

		ippsSet_32f(xdf->valueBefore, xd->matd1, nxz);
		ippsSet_32f(xff->valueBefore, xd->matf1, nxz);

		if(xdf->typeProfile == FIELD_PROFILE_CONSTANT){
			ippsSet_32f(xdf->valueBefore, xd->matd2, nxz);
		}else{
			ippsSet_32f(xdf->valueAfter, xd->matd2, nxz);
		}

		if(xff->typeProfile == FIELD_PROFILE_CONSTANT){
			ippsSet_32f(xff->valueBefore, xd->matf2, nxz);
		}else{
			ippsSet_32f(xff->valueAfter, xd->matf2, nxz);
		}
	}

	if(xdf->typeProfile == FIELD_PROFILE_PROFILE){
		if(aD->ny != aD->profd_nx){
			sprintf(aD->message,"height of a slice (%u) should =  width of dark profile (%u)", aD->ny, aD->profd_nx);
			printError(aD);
			return false;
		}
	}

	xd->veci = ippsMalloc_32f(aD->nt);
	if(xd->veci == NULL){
		printError(aD,"can not allocate memory for input data");
		return false;
	}
	if(xid->restrictions == INPUT_RESTRICTIONS_YES){
		if(aD->nb == 32){
			printInfo(aD,"restrictions will not be applied");
		}else{
			vmin = xid->valueMin;
			vmax = xid->valueMax;
			sprintf(aD->message, "Minimal (%f) and maximal (%f) values", vmin, vmax);
			printInfo(aD);
			if(vmax - vmin < 0.00001f){
				printError(aD, "wrong values");
				return false;
			}
		}
	}
	xd->vecd = ippsMalloc_32f(aD->nt);
	xd->vecf = ippsMalloc_32f(aD->nt);
	if(xd->vecd == NULL || xd->vecf == NULL){
		printError(aD,"can not allocate memory for dark/flat slices");
		return false;
	}

	if(xdf->typeProfile == FIELD_PROFILE_PROFILE){
		if(aD->profd_ny <= aD->row_last){
			sprintf(aD->message, "wrong height of the dark profile (%u), should be greater than the last slice (%u)", aD->profd_ny, aD->row_last);
			printError(aD);
			return false;
		}
		if(aD->profd_nx != aD->ny){
			sprintf(aD->message, "wrong width of the dark profile (%u), should  = the height of the slice (%u)", aD->profd_nx, aD->ny);
			printError(aD);
			return false;
		}
	}

	if(xff->typeProfile == FIELD_PROFILE_PROFILE){
		if(aD->proff_ny <= aD->row_last){
			sprintf(aD->message, "wrong height of the flat profile (%u), should be greater than the last slice (%u)", aD->proff_ny, aD->row_last);
			printError(aD);
			return false;
		}
		if(aD->proff_nx != aD->ny){
			sprintf(aD->message, "wrong width of the flat profile (%u), should  = the height of the slice (%u)", aD->proff_nx, aD->ny);
			printError(aD);
			return false;
		}
	}

	xd->vect1 = ippsMalloc_32f(aD->nx);
	xd->vect2 = ippsMalloc_32f(aD->nx);
	if(xd->vect1 == NULL || xd->vect2 == NULL){
		printError(aD,"can not allocate memory (temp vectors)");
		return false;
	}
	return true;
}



bool read_one_slice(allData *aD){
	long fs;
	float vmin, vmax, sd;

	xtiff *ti;
	FILE *ff;
	xData *xd;
	xInputData *xid;

	xd = aD->data;
	xid = aD->hms->fbp->inputData;

	sprintf(aD->ifile,aD->itemplate,aD->row);
	printTag(aD,"InputFile",aD->ifile);
	
	ff = fopen(aD->ifile,"rb");
	if(ff == NULL){
		printError(aD,"can not open the input file");
		return false;
	}
	
	fseek(ff,0,SEEK_END);
	fs = ftell(ff);

	printTag(aD,"FileSize",fs,"file size (bytes)");
	if(fs < xid->memorySizeMin || fs > xid->memorySizeMax){
		printError(aD,"wrong file size");
		return false;
	}
	fseek(ff,0,SEEK_SET);
	fread(xd->vecio,sizeof(Ipp8u),fs,ff);
	fclose(ff); ff = NULL;

	ti = new xtiff();
	if(!ti->getInfo(xd->vecio,fs,&(aD->temp_nx),&(aD->temp_ny),&(aD->temp_nb))){
		printError(aD,"can not read the TIFF file");
		return false;
	}
	printTag(aD,"FileWidth",aD->temp_nx);
	printTag(aD,"FileHeight",aD->temp_ny);
	printTag(aD,"FileBits",aD->temp_nb);
	
	if(aD->isFirstSlice){
		if(!read_one_slice_fs(aD))return false;
	}

	if(aD->temp_nx != aD->nx){
		printError(aD,"slices have different widths");
		return false;
	}
	if(aD->temp_ny != aD->ny){
		printError(aD,"slices have different heights");
		return false;
	}
	if(aD->temp_nb != aD->nb){
		printError(aD,"slices have different number of bits");
		return false;
	}
	
	vmin = xid->valueMin;
	vmax = xid->valueMax;

	if(!ti->read_float(xd->vecio,fs,xd->veci)){
		sprintf(aD->message,"can not read TIFF file (%s)",ti->err_message);
		delete ti;
		printError(aD);
		return false;
	}

	delete ti;
	
	if(xid->restrictions == INPUT_RESTRICTIONS_YES && aD->nb != 32){
		sd = vmax - vmin;
		if(aD->nb == 8){
			sd /= 255.0f;
		}else{
			sd /= 65535.0f;
		}
		ippsMulC_32f_I(sd,xd->veci,aD->nt);
		ippsAddC_32f_I(vmin,xd->veci,aD->nt);
	}

	return true;
}



bool print_output_data(allData *aD){
	int nx, ny, nb, nt;
	long ofile_size;
	FILE *ff;
	xtiff *to;
	Ipp32f minv, maxv, dv, mv;
	xData *xd;

	xd = aD->data;

	printTagStart(aD,"PrintOutputData");

	sprintf(aD->ofile,aD->otemplate,aD->row);
	printTag(aD,"File",aD->ofile);
	
	ff = fopen(aD->ofile,"wb");
	if(ff == NULL){
		printError(aD,"can not open the output file");
		return false;
	}

	if(aD->hms->fbp->outputData->bitsType == OUTPUT_BITS_GIVEN){
		nb = aD->hms->fbp->outputData->bits;
	}else{
		nb = aD->nb;
	}
	nx = aD->nx_out;
	ny = aD->ny_out;
	nt = nx * ny;

/*	for(int j=1; j<ny; j++){
		ippsSub_32f_I(xd->veco, xd->veco +j*nx,nx);
	}
	ippsZero_32f(xd->veco,nx);
*/
	if(aD->hms->fbp->outputData->state == OUTPUT_STATE_INTENSITY){
		ippsMulC_32f_I(-1.0,xd->veco,nt);
		ippsExp_32f_I(xd->veco,nt);
	}

	to = new xtiff();
	ofile_size = to->get_file_size(nx,ny,nb);
	printTag(aD,"FileSize",ofile_size,"file size (bytes)");
	if(ofile_size >aD->hms->fbp->inputData->memorySizeMax){
		printError(aD,"wrong file size");
		return false;
	}
	
	if(nb < 32){
		if(aD->hms->fbp->outputData->restrictions == OUTPUT_RESTRICTIONS_YES){
			minv = aD->hms->fbp->outputData->valueMin;
			maxv = aD->hms->fbp->outputData->valueMax;
			if(minv > maxv){
				printTag(aD,"MinimalValue",minv);
				printTag(aD,"MaximalValue",maxv);
				printError(aD,"wrong values");
				return false;
			}
		}else{
			ippsMinMax_32f(xd->veco,nt,&minv,&maxv);
		}
		dv = maxv - minv;
		if(dv < 1.e-6){
			mv = 0.0f;
		}else{
			if(nb == 8){
				mv = 255.0f/dv;
			}else{
				mv = 65535.0f/dv;
			}
		}
		printTag(aD,"MinimalValue",minv);
		printTag(aD,"MaximalValue",maxv);
		ippsSubC_32f_I(minv,xd->veco,nt);
		ippsMulC_32f_I(mv,xd->veco,nt);
	}

	if(!to->print_float(xd->vecio,xd->veco,nx,ny,nb)){
		sprintf(aD->message,"output TIFF (%s)",to->err_message);
		printError(aD);
		return false;
	}

	fwrite(xd->vecio,sizeof(Ipp8u),ofile_size,ff);
	fclose(ff); ff = NULL;

	delete to;
	printTagEnd(aD);
	return true;
}


