#include "defs.h"

xXrm2tif::xXrm2tif(){
	width = 0;
	height = 0;
	NOI = 1;
	NOD = 5;
	bitsPerData = -1;
	bitsPerRef = -1;
	pixelSize = 0.0;
	d2a = 0.0;
	s2a = 0.0;
}

bool main_xrm2tif(allData *aD){
	char tname[MAX_FOLDER], finfo[MAX_FOLDER], fraw[MAX_FOLDER];
	char iname[MAX_FOLDER];
	xXrm2tif* x2t; 
	ole2cd *cd;

	x2t = aD->hms->x2t; 
	aD->fp = fopen(x2t->logFile,"w");
	if(aD->fp == NULL){
		printf("Can not open log file \"%s\".\n",x2t->logFile);
		return false;
	}
	printTagStart(aD,"logXrm2tif");
	print_global_time_stamp(aD);
	
	sprintf(tname,"%s",x2t->folder);

	if(!ifFolderExists(aD,tname)){
		if(!createFolderCycle(aD,tname)){
			return false;
		}
	}
	
	sprintf(finfo,"%s/info",tname);
	sprintf(fraw,"%s/raw",tname);
	if(!createFolderCycle(aD,finfo))return false;
	if(!createFolderCycle(aD,fraw))return false;

	cd = new ole2cd();
	sprintf(iname,"%s",x2t->inputFile);
	cd->ff = fopen(iname,"rb");
	if(cd->ff == NULL){
		fprintf(aD->fp,"Error: can not open input file \"%s\".\n", iname);
		return false;
	}

	sprintf(aD->message,"Input file \"%s\" has been opened",iname);
	printInfo(aD);

	if(!cd->read_header()){delete cd; return false;}
	if(!cd->sort_dir()){delete cd; return false;}
	if(!cd->read_types(aD)){delete cd; return false;}
	if(!cd->print_xml(aD,finfo)){delete cd; return false;}
	if(!cd->fill_MainParam(aD)){delete cd; return false;}
	if(!cd->print_MainParam(aD, finfo)){delete cd; return false;}

	if(x2t->output_type == XRM_OUTPUT_ALL){
		if(x2t->bitsPerRef>0){
			if(!cd->print_ref_data(aD)){delete cd; return false;}
		}else{
			printInfo(aD,"There is no reference file");
		}

		if(x2t->bitsPerData>0){
			if(!cd->print_data(aD,fraw)){delete cd; return false;}
		}else{
			printInfo(aD,"There are no images");
		}
	}
	delete cd;
	
	return true;
}

