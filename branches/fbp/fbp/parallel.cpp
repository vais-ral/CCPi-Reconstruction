#include "defs.h"
#include "time_stamp.h"


bool main_fbp(allData *aD){
	xFBP* xfbp; 
	xData *xd;
	xOutputData *xod;
	xInputData *xid;

	xd = aD->data;
	xid = aD->hms->fbp->inputData;
	xod = aD->hms->fbp->outputData;
	xfbp = aD->hms->fbp; 

	timestamp("Entering main_fbp",4);


	aD->fp = fopen(xfbp->logFile,"w");
	if(aD->fp == NULL){
		printf("Can not open log file \"%s\".\n",xfbp->logFile);
		return false;
	}
	printTagStart(aD,"logFBP");
	print_global_time_stamp(aD);

	printSettingsXml_FBP(aD);
	
	if(!allocate_memory(aD))return false;
	if(!read_fdf_data(aD))return false;
	if(!check_rows(aD))return false;
	if(!form_files_names(aD))return false;

	for(unsigned int row = aD->row_first; row<= aD->row_last; row+= aD->row_step){
		if(row>aD->row_first){
			aD->isFirstSlice = false;
		}
		sprintf(aD->message,"Slice %04u",row);
		timestamp(aD->message,4);
		printTagStart(aD,aD->message);
		if(row > aD->row_first) aD->isFirstSlice = false;
		aD->row = row;
		timestamp("calling read_one_slice",4);
		if(!read_one_slice(aD))return false;
		timestamp("returned from read_one_slice",4);

		if(xid->type == INPUT_TYPE_INTENSITY){
			if(!ff_correction(aD))return false;
		}
		timestamp("returned from ff_correction",4);

		if(xod->type == OUTPUT_TYPE_FFCORRECTION){
			xd->veco = xd->veci;
			aD->nx_out = aD->nx;
			aD->ny_out = aD->ny;
			print_output_data(aD);
			printTagEnd(aD);
			continue;
		}
		

		
		timestamp("Calling preprocessing",4);
		printTagStart(aD,"Preprocessing");
		if(!high_peaks_before(aD))return false;
		timestamp("returned from high_peaks_before",4);
		if(!ring_artefacts_removal(aD))return false;
		timestamp("returned from ring_artefacts_removal",4);
		if(!intensity_norm(aD))return false;
		timestamp("returned from intensity_norm",4);
		printTagEnd(aD);

		if(xod->type == OUTPUT_TYPE_PREPROCESSED){
			xd->veco = xd->veci;
			aD->nx_out = aD->nx;
			aD->ny_out = aD->ny;
			print_output_data(aD);
			printTagEnd(aD);
			continue;
		}
		timestamp("Returned from preprocessing",4);
		timestamp("Calling transform_sinogram",4);
		
		if(!transform_sinogram(aD))return false;

		timestamp("Returned from transform_sinogram",4);

		if(xod->type == OUTPUT_TYPE_TRANSFORMED){
			xd->veco = xd->vta;
			aD->nx_out = aD->ta_nx;
			aD->ny_out = aD->ta_ny;
			print_output_data(aD);
			printTagEnd(aD);
			continue;
		}
		
		timestamp("Calling fbp_gpu",4);
		if(!fbp_gpu(aD))return false;
		timestamp("Returned from fbp_gpu",4);
	/*	
		aD->nx_out = aD->gi->pol_a;
		aD->ny_out = aD->gi->pol_r;
		xd->veco = xd->vecPol;
*/
		
		aD->nx_out = aD->gi->mxo;
		aD->ny_out = aD->gi->myo;
		xd->veco = xd->vto;
		
		
		timestamp("calling print_output_data",4);
		print_output_data(aD);
		timestamp("Returned from print_output_data",4);
		printTagEnd(aD);
	}

	timestamp("Leaving main_fbp",4);

	return true;
}
