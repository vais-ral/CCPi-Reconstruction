/* ----------------------------------------------------------------------
 * allFDK.cpp()
 * ----------------------------------------------------------------------
 * FDK reconstruction routines for HMXIF code. Calls the kernel launching
 * code in the .cu file.
 * * 
 * ----------------------------------------------------------------------
 */

#include "defs.h"


bool read_vector_xml(allData *aD, char *inputFile, Ipp32f *vec){
	char tagName[100], childTag[100], value[MAX_FOLDER];
	int pos1, pos2, pos;
	int s1, s2, s3, s4;
	int ch1, ch2, ch3, ch4;
	int n, nm, u2, u3, lenf;
	bool t, isEmptyTag;	
	Ipp8u str_comm_left[5], str_comm_right[5], str_xml_left[3], str_xml_right[3], str_tag_left[2], str_tag_right[2], str_space[2];
	FILE *ff;
	Ipp8u *vt, *vf;

	sprintf((char *)str_comm_left,"<!--");
	sprintf((char *)str_comm_right,"-->");
	sprintf((char *)str_tag_left,"<");
	sprintf((char *)str_tag_right,">");
	sprintf((char *)str_space," ");
	sprintf((char *)str_xml_left,"<?");
	sprintf((char *)str_xml_right,"?>");

	printTag(aD,"InputFile",inputFile);

	ff = fopen(inputFile,"r");
	if(ff == NULL){
		printError(aD,"can not open the file");
		return false;
	}
	fseek(ff,0,SEEK_END);
	lenf = ftell(ff);
	fseek(ff,0,SEEK_SET);
	printTag(aD,"FileSize",lenf,"file size in bytes");

	vt = ippsMalloc_8u(lenf);
	vf = ippsMalloc_8u(lenf);
	if(vf == NULL || vt == NULL){
		printError(aD,"can not allocate memory to read the file");
		return false;
	}

	fread(vf,sizeof(Ipp8u),lenf,ff);
	fclose(ff); ff = NULL;
	
	ippsAddC_8u_ISfs(128, vf, lenf, 0);
	ippsSubC_8u_ISfs(128, vf, lenf, 0);
	ippsSubC_8u_ISfs(32, vf, lenf, 0);
	ippsAddC_8u_ISfs(32, vf, lenf, 0);

	ippsReplaceC_8u(vf, vt, lenf, 127, 32);
	ippsCopy_8u(vt,vf,lenf);
	ippsFree(vt); vt = NULL;
		
	for(;;){
		ippsFind_8u(vf,lenf,str_xml_left, 2, &pos1);
		if(pos1 <0)break;
		ippsFind_8u(vf,lenf,str_xml_right, 2, &pos2);
		if(pos2<pos1){
			printError(aD,"xml file, header");
			return false;
		}
		ippsRemove_8u_I(vf, &lenf, pos1, pos2-pos1+2);
	}

	for(;;){
		ippsFind_8u(vf,lenf,str_comm_left, 4, &pos1);
		if(pos1 <0)break;
		ippsFind_8u(vf,lenf,str_comm_right, 3, &pos2);
		if(pos2<pos1){
			printError(aD,"xml file, comments");
			return false;
		}
		ippsRemove_8u_I(vf, &lenf, pos1, pos2-pos1+3);
	}

	pos = 0;
	for(;;){
		ippsFind_8u(vf+pos,lenf-pos,str_tag_left, 1, &pos2);
		if(pos2 <0)break;
		pos+=pos2;
		ippsFind_8u(vf+pos,lenf-pos,str_tag_right, 1, &pos1);
		ippsFind_8u(vf+pos,lenf-pos,str_space, 1, &pos2);
		if(pos2 <0 || pos2 >pos1){
			pos+=pos1;
			continue;
		}
		ippsRemove_8u_I(vf, &lenf, pos+pos2, pos1-pos2);
		pos++;
	}

	if(vf[0] == 32){
		ippsFind_8u(vf,lenf,str_tag_left, 1, &pos1);
		ippsRemove_8u_I(vf, &lenf, 0, pos1);
	}
	if(vf[lenf-1] == 32){
		ippsFindRev_8u(vf,lenf,str_tag_right, 1, &pos1);
		ippsRemove_8u_I(vf, &lenf, pos1+1, lenf-pos1-1);
	}

	pos = 0;
	do{
		ippsFind_8u(vf+pos,lenf,str_tag_right, 1, &pos1);
		if(pos1<0)break;
		pos+=pos1;
		pos2 = pos;
		for(int i=pos+1;i<lenf;i++){
			if(vf[i] != 32) break;
			pos2++;
		}
		if(pos2 > pos)ippsRemove_8u_I(vf, &lenf, pos+1, pos2-pos);

		
		pos++;
	}while(pos<lenf);

	pos = lenf-1;
	
	do{
		ippsFindRev_8u(vf,pos,str_tag_left, 1, &pos1);
		if(pos1<0)break;
		pos=pos1;
		pos2 = pos;
		for(int i=pos-1;i>=0;i--){
			if(vf[i] != 32) break;
			pos2--;
		}
		if(pos2 < pos){
			ippsRemove_8u_I(vf, &lenf, pos2, pos-pos2);
		}
		pos--;
	}while(pos>=0);
	
	t = findNextTag(vf, 0, lenf-1, &s1, &s2, tagName, &isEmptyTag);
	if(!t){
		printError(aD,"can not find the root start-tag");
		return false;
	}
	
	t = findEndTag(vf, s2, lenf-1, &s3, &s4, tagName);
	if(!t){
		printError(aD,"can not find the root end-tag");
		return false;
	}

	
	u2 = s2;
	u3 = s3;
	
	for(;;){
		t = findNextTag(vf, s2, s3, &ch1, &ch2, childTag, &isEmptyTag);		
		if(!t)break;
		printInfo(aD,childTag);
		t = findEndTag(vf, ch2, s3, &ch3, &ch4, childTag);
		if(!t){
			sprintf(aD->message,"can not find the child end-tag \"%s\"",childTag);
			printError(aD);
			return false;
		}
		u2 = s2;
		u3 = s3;
		s2 = ch2 + 1;
		s3 = ch3 - 1;
	}
	s2 = u2;
	s3 = u3;

	nm = aD->hms->af->inputImageLast + 1;
	
	n = 0;
	printTagStart2(aD,"Vector");
	for(;;){
		t = findNextTag(vf, s2, s3, &ch1, &ch2, childTag, &isEmptyTag);
		if(!t)break;
		t = findEndTag(vf, ch2, s3, &ch3, &ch4, childTag);
		if(!t){
			sprintf(aD->message,"can not find the child end-tag \"%s\"",childTag);
			printError(aD);
			return false;
		}
		
		getValue(vf,ch2,ch3,value);
		vec[n] = float(atof(value));
		sprintf(aD->message,"%i",n);
		printTag(aD,"Value",vec[n],aD->message);
		s2 = ch4+1;
		n++;
		if(n > nm) break;	
	}
	printTagEnd2(aD);

	for(int i = n; i < nm; i++){
		vec[i] = vec[n-1];
	}
	ippsFree(vt); vt = NULL;
	ippsFree(vf); vf = NULL;
	return true;

}

bool motion_correction_af(allData *aD){
	xAllFDK* af;
	xData *xd;
	int nm;
	
	printTagStart2(aD,"MotionCorrection");
	xd = aD->data;
	af = aD->hms->af;
	nm = af->inputImageLast+1;

	printTagStart2(aD,"HorizontalShifts");

	xd->shift_x = ippsMalloc_32f(nm);
	if(xd->shift_x == NULL){
		printError(aD,"can not allocate memory (horizontal shifts)");
		return false;
	}
	
	if(af->xShiftsType == SHIFTS_NO){
		printInfo(aD,"No shifts");
		ippsZero_32f(xd->shift_x,nm);
	}else{
		if(!read_vector_xml(aD,af->xShiftsFile,xd->shift_x))return false;
	}
	printTagEnd2(aD);

	printTagStart2(aD,"VerticalShifts");

	xd->shift_y = ippsMalloc_32f(nm);
	if(xd->shift_y == NULL){
		printError(aD,"can not allocate memory (vertical shifts)");
		return false;
	}
	
	if(af->xShiftsType == SHIFTS_NO){
		printInfo(aD,"No shifts");
		ippsZero_32f(xd->shift_y,nm);
	}else{
		if(!read_vector_xml(aD,af->yShiftsFile,xd->shift_y))return false;
	}
	printTagEnd2(aD);
	
	printTagEnd2(aD);
	return true;
}


bool ring_artefacts_af(allData *aD){
	int nx, ny, nz, nt;
	int ii, jj;
	int t;
	Ipp64f mdiv, alpha, tau, paramN, paramR;
	xAllFDK* af;
	xData *xd;
	
	printTagStart(aD,"RingArtefacts");
	af = aD->hms->af;
	xd = aD->data;
	
	nx = aD->nx;
	ny = aD->ny;
	nz = aD->nz;
	nt= nx*ny;
	
	paramN = double(af->ringArtefactsParamN);
	paramR = double(af->ringArtefactsParamR);
	alpha = (1+paramN)*paramR;

	fprintf(aD->fp,"Parameter 1 (ratio of norms): %f.\n",paramN);
	fprintf(aD->fp,"Parameter 2 (regularisation): %f.\n",paramR);
	if(paramN <-1e-10 || paramR <1e-8){
		fprintf(aD->fp,"Error: wrong parameters.\n");
		return false;
	}

	t = af->ringArtefactsType;
	
	if(t == RING_ARTEFACTS_COLUMN || t == RING_ARTEFACTS_AML){
		xd->vec32 = ippsMalloc_32f(nx);
		xd->vec64 = ippsMalloc_64f(nx);
		xd->vec_res = ippsMalloc_64f(nx);
		if(xd->vec32 == NULL || xd->vec64 == NULL || xd->vec_res == NULL){
			fprintf(aD->fp,"Error: can not allocate memory (ring artefacts).\n");
			return false;
		}
	}

	if(t == RING_ARTEFACTS_AML){
		Ipp64f *vec1, *vec2;
		vec1 = ippsMalloc_64f(1);
		vec2 = ippsMalloc_64f(1);
		Ipp64f *vec_exp = ippsMalloc_64f(2*nx+2);
		xd->mata = ippsMalloc_64f(nx*nx);
		if(vec1 == NULL || vec2 == NULL || xd->mata == NULL || vec_exp == NULL){
			fprintf(aD->fp,"Error: can not allocate memory (ring artefacts, matrix).\n");
			return false;
		}
		
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
				xd->mata[i*nx+j] = vec_exp[abs(jj-ii)]* (1.0+vec_exp[2*nx-2*jj+1]) * (1.0+vec_exp[2*ii-1]);
			}
		}
		ippsDivC_64f_I(mdiv*(1.0-vec_exp[2*nx]),xd->mata,nx*nx);

		ippsFree(vec1); vec1 = NULL;
		ippsFree(vec2); vec2 = NULL;
		ippsFree(vec_exp); vec_exp = NULL;
	}
	

	switch(t){
		case RING_ARTEFACTS_NO:
			fprintf(aD->fp,"No ring artefact suppression.\n");
			break;
		case RING_ARTEFACTS_COLUMN:
			fprintf(aD->fp,"A mean value is subtracted from each column.\n");
			for(int k=0;k<ny;k++){
				ippsZero_64f(xd->vec_res,nx);
				for(int i=0;i<nz;i++){
					ippsConvert_32f64f(xd->veci+i*nt+k*nx,xd->vec64,nx);
					ippsAdd_64f_I(xd->vec64,xd->vec_res,nx);
				}
				ippsDivC_64f_I(double(nz),xd->vec_res,nx);
				ippsConvert_64f32f(xd->vec_res,xd->vec32,nx);
				for(int i=0;i<nz;i++){
					ippsSub_32f_I(xd->vec32,xd->veci+i*nt+k*nx,nx);
				}
			}
			break;
		case RING_ARTEFACTS_AML:
			fprintf(aD->fp,"Suppression based on paper \"Applied Mathematics Letters, v. 23, p. 1489\".\n");
			for(int k=0;k<ny;k++){
				ippsZero_64f(xd->vec_res,nx);
				for(int i=0;i<nz;i++){
					ippsConvert_32f64f(xd->veci+i*nt+k*nx,xd->vec64,nx);
					ippsAdd_64f_I(xd->vec64,xd->vec_res,nx);
				}
				ippsDivC_64f_I(double(nz),xd->vec_res,nx);

				ippsZero_64f(xd->vec64,nx);
				ippsMulC_64f(xd->vec_res,paramN*paramR,xd->vec64,nx);
				ippsAdd_64f_I(xd->vec_res,xd->vec64,nx-1);
				ippsSub_64f_I(xd->vec_res+1,xd->vec64,nx-1);
				ippsAdd_64f_I(xd->vec_res+1,xd->vec64+1,nx-1);
				ippsSub_64f_I(xd->vec_res,xd->vec64+1,nx-1);

				for(int i=0;i<nx;i++){
					ippsDotProd_64f(xd->vec64,xd->mata+i*nx,nx,xd->vec_res + i);
				}
				ippsConvert_64f32f(xd->vec_res,xd->vec32,nx);
			
				for(int i=0;i<nz;i++){
					ippsSub_32f_I(xd->vec32,xd->veci+i*nt+k*nx,nx);
				}
			}
			break;
		default:
			break;
	}

	printTagEnd(aD);
	
	return true;
}


bool ff_correction_af(allData *aD){
	printf("Flat field correction...\n");
	xData *xd = aD->data;
	int nx, ny, nz, nt;

	nx = aD->nx;
	ny = aD->ny;
	nz = aD->nz;
	nt= nx*ny;

	ippsSub_32f_I(xd->matd1,xd->matf1,nt);
	ippsThreshold_LT_32f_I(xd->matf1,nt,0.1f);
	for(int i=0;i<nz;i++){
		ippsSub_32f_I(xd->matd1,xd->veci+i*nt,nt);
	}
	ippsThreshold_LT_32f_I(xd->veci,nz*nt,0.1f);
	for(int i=0;i<nz;i++){
		ippsDiv_32f_I(xd->matf1,xd->veci+i*nt,nt);
	}
	ippsLn_32f_I(xd->veci,nz*nt);
	ippsMulC_32f_I(-1.0,xd->veci,nz*nt);
	
	printf("Done\n\n");
	return true;
}


bool read_input_data_af(allData *aD){
	char ftemplate[MAX_FOLDER], fname[MAX_FOLDER];
	unsigned int nx, ny, nz, nb;
	unsigned int ipf, nh, np;
	int old_file, new_file;
	unsigned int imf, iml, ims;
	int nnz;
	unsigned int mx, my, mb;
	float v;
	xAllFDK* af;
	xData *xd;
	
	printTagStart(aD,"InputData");
	af = aD->hms->af;
	xd = aD->data;
	
	if(af->inputFileFirst <0){
		sprintf(aD->message,"wrong input file number (%i)",af->inputFileFirst);
		printError(aD);
		return false;
	}
	if(af->inputFileStep <1){
		sprintf(aD->message,"wrong step for input files (%i)",af->inputFileStep);
		printError(aD);
		return false;
	}
	sprintf(ftemplate,"%s/%s%%0%ii%s.%s",af->inputFolder,af->inputPrefix,af->inputNOD,af->inputSuffix,af->inputExtension);

	sprintf(fname,ftemplate,af->inputFileFirst);
	gtiff *ti = new gtiff();
	if(!ti->getInfo(fname,&nx,&ny,&nb)){
		sprintf(aD->message,"can not get information about file \"%s\"",fname);
		printError(aD);
		return false;
	}
	printTagStart2(aD,"FilesInfo");
	printTag(aD,"Width",nx);
	printTag(aD,"Height",ny);
	printTag(aD,"NumberOfBits",nb);
	printTagEnd2(aD);
	//fprintf(aD->fp,"Input files: width (%i), height (%i), number of bits (%i).\n",nx,ny,nb);
	delete ti;

	ipf = af->inputImagesPerFile;

	printTag(aD,"ImagesPerFile",ipf);
	if(ipf < 1 || ipf > ny){
		printError(aD,"wrong number");
		return false;
	}

	if(ny%ipf != 0){
		fprintf(aD->fp,"Error: wrong number.\n");
		return false;
	}

	ny /= ipf;

	printTag(aD,"HeightOfEachImage",ny);
	
	imf = af->inputImageFirst;
	iml = af->inputImageLast;
	ims = af->inputImageStep;
	printTag(aD,"ImageFirst",imf);
	printTag(aD,"ImageLast",iml);
	printTag(aD,"ImageStep",ims);
	if(ims < 1){
		printError(aD,"wrong number");
		return false;
	}
	if(iml< imf){
		printError(aD,"wrong numbers");
		return false;
	}

	nz = (iml-imf)/ims+1;
	sprintf(aD->message,"Total number of images to be read: %u",nz);
	printInfo(aD);
	
	old_file = -1;

	xd->veci = ippsMalloc_32f(nx*ny*nz);
	if(xd->veci == NULL){
		sprintf(aD->message,"can not allocate memory to read input files (%u x %u x %u)",nx,ny,nz);
		printError(aD);
		return false;
	}
	xd->vect1 = ippsMalloc_32f(nx*ny*ipf);
	if(xd->vect1 == NULL){
		printError(aD,"can not allocate memory to read input files");
		return false;
	}

	ti = new gtiff();
	
	nnz = 0;
	printTagStart2(aD,"ReadInputFiles");
	for(unsigned int i = imf; i <= iml; i += ims){
		new_file = af->inputFileFirst + af->inputFileStep *((i-imf)/ipf);
		if(new_file != old_file){
			sprintf(fname,ftemplate,new_file);
			printTag(aD,"InputFile",fname);
			if(!ti->read_float(fname,xd->vect1)){
				printError(aD,"can not read data from the file");
				return false;
			}
			nh=ti->ImageLength/ny;
			if(ti->ImageLength%ny != 0 || nh == 0){
				sprintf(aD->message,"wrong height of the image (%u)",ti->ImageLength);
				printError(aD);
				return false;
			}
			old_file = new_file;
		}
		np = i%ipf;
		if(np+1> nh){
			sprintf(aD->message,"wrong number of chunks (%u) in the file",nh);
			printError(aD);
			return false;
		}
		ippsCopy_32f(xd->vect1+np*nx*ny,xd->veci+nnz*nx*ny,nx*ny);
		nnz++;
	}
	printTagEnd2(aD);
	delete ti;

	
	if(af->inputRestrictions == INPUT_RESTRICTIONS_YES && nb < 32){
		fprintf(aD->fp,"Input restrictions: %f (min), %f (max).\n",af->inputMin,af->inputMax);
		if(af->inputMin > af->inputMax - 1.e-6){
			fprintf(aD->fp,"Error: wrong values.\n");
			return false;
		}
		v = af->inputMax - af->inputMin;
		if(nb == 8){
			v/= 255.0;
		}else{
			v/=65535.0;
		}
		ippsMulC_32f_I(v,xd->veci,nx*ny*nz);
		ippsAddC_32f_I(af->inputMin,xd->veci,nx*ny*nz);
	}

	xd->matd1 = ippsMalloc_32f(nx*ny);
	xd->matf1 = ippsMalloc_32f(nx*ny);
	if(xd->matd1 == NULL || xd->matf1 == NULL){
		printError(aD,"can not allocate memory to read flat/dark fields");
		return false;
	}

	ti = new gtiff();
	if(af->darkType == FIELD_TYPE_VALUE){
		fprintf(aD->fp,"Dark field, user value (%f).\n",af->darkValue);
		if(af->darkValue < -1.e-6){
			printError(aD,"wrong value");
			return false;
		}
		ippsSet_32f(af->darkValue,xd->matd1,nx*ny);
	}else{
		fprintf(aD->fp,"Dark field file: \"%s\".\n",af->darkName);
		if(!ti->getInfo(af->darkName,&mx,&my,&mb)){
			fprintf(aD->fp,"Error: can not get information about the file \"%s\".\n",af->darkName);
			return false;
		}
		fprintf(aD->fp,"Dark field file: width (%u), height (%u), number of bits (%u).\n",mx,my,mb);
		if(nx != mx || ny != my){
			printError(aD,"wrong size");
			return false;
		}
		if(!ti->read_float(af->darkName,xd->matd1)){
			printError(aD,"can not read data from the dark field file");
			return false;
		}
		if(mb == 8 || mb == 16){
			fprintf(aD->fp,"Dark field: min (%f) and max (%f) values.\n",af->darkMin,af->darkMax);
		
			if(af->darkMin > af->darkMax - 1.e-6){
				printError(aD,"wrong values");
				return false;
			}
			v = af->darkMax - af->darkMin;
			if(mb == 8){
				v/= 255.0f;
			}else{
				v/=65535.0f;
			}
			ippsMulC_32f_I(v,xd->matd1,nx*ny);
			ippsAddC_32f_I(af->darkMin,xd->matd1,nx*ny);
		}
	}
	delete ti;

	ti = new gtiff();
	if(af->flatType == FIELD_TYPE_VALUE){
		fprintf(aD->fp,"Flat field, user value (%f).\n",af->flatValue);
		if(af->flatValue < 1.e-3){
			printError(aD,"wrong value");
			return false;
		}
		ippsSet_32f(af->flatValue,xd->matf1,nx*ny);
	}else{
		fprintf(aD->fp,"Flat field file: \"%s\".\n",af->flatName);
		if(!ti->getInfo(af->flatName,&mx,&my,&mb)){
			fprintf(aD->fp,"Error: can not get information about the file \"%s\".\n",af->flatName);
			return false;
		}
		fprintf(aD->fp,"Flat field file: width (%u), height (%u), number of bits (%u).\n",mx,my,mb);
		if(nx != mx || ny != my){
			printError(aD,"wrong size");
			return false;
		}
		if(!ti->read_float(af->flatName,xd->matf1)){
			printError(aD,"can not read data from the flat field file");
			return false;
		}
		if(mb == 8 || mb == 16){
			fprintf(aD->fp,"Flat field: min (%f) and max (%f) values.\n",af->flatMin,af->flatMax);
		
			if(af->flatMin > af->flatMax - 1.e-6){
				printError(aD,"wrong values");
				return false;
			}
			v = af->flatMax - af->flatMin;
			if(mb == 8){
				v/= 255.0f;
			}else{
				v/=65535.0f;
			}
			ippsMulC_32f_I(v,xd->matf1,nx*ny);
			ippsAddC_32f_I(af->flatMin,xd->matf1,nx*ny);
		}
	}
	delete ti;

	aD->nx = nx;
	aD->ny = ny;
	aD->nz = nz;
	
	printTagEnd(aD);
	return true;
}

bool main_allfdk(allData *aD){
//	int numProj;
	int mw, ew;
	int outfile;
	int bs, ext_left, ext_right, ext_num, ext_rest;
	int nx, ny, nwx, order, nz;
	int nbo;
	int outb;
	IppsFFTSpec_R_32f* spec;
	Ipp32f *vt;
	float angular_step;
	float slf, sll, sls;
	char tname[MAX_FOLDER];
	char otemplate[MAX_FOLDER], ofile[MAX_FOLDER];
	xAllFDK* af; 
	PolarToCart *PC;
	MatPar *matp;
	OtherParam *other;
	Ipp32f sa;
	Ipp8u *buf;
	gtiff *to;




	af = aD->hms->af; 
	aD->fp = fopen(af->logFile,"w");
	if(aD->fp == NULL){
		printf("Can not open log file \"%s\".\n",af->logFile);
		return false;

	}
	printTagStart(aD,"logAllFDK");
	print_global_time_stamp(aD);

	if(!read_input_data_af(aD))return false;
	if(!ff_correction_af(aD))return false;
	if(!ring_artefacts_af(aD))return false;

	if(!motion_correction_af(aD))return false;

	sprintf(tname,"%s",af->outputFolder);

	if(!ifFolderExists(aD,tname)){
		if(!createFolderCycle(aD,tname)){
			return false;
		}
	}

//	numProj = aD->nz;
	angular_step = af->angleStep;
	
	printTag(aD,"AngularStep", angular_step);
	if(angular_step <0.01){
		printError(aD,"wrong number");
		return false;
	}

	other = new OtherParam();
	other->w_a = 1; other->w_r = 4; other->w_s = 1;
	other->images_to_process = int(180.0/af->angleStep);
	other->num_images = 2*other->images_to_process;

	PC = new PolarToCart();
	if(!fillPC(PC,aD,other)){
		printError(aD,"can not form PC structure");
		return false;
	}

	if(!PC->formCudaMatrix(aD)){
		printError(aD,"can not allocate memory (polar)");
		return false;
	}

	printTag(aD,"Rmin",PC->rmin);
	printTag(aD,"Rmax",PC->rmax);
	
	mw = aD->nx-1;
	ew = 1;
	while(mw>1){
		ew++;
		mw/=2;
	}
	other->eWidth = 1 << ew;
	formRowFilter(aD,other);
	
	nz = aD->nz;
	nx = aD->nx;
	nwx = other->eWidth;
	order = ew;
	ext_num = __min(__max(5,10),200);
	ext_rest=nwx-nx;
	ext_left=(ext_rest/2);
	ext_right=ext_rest - ext_left;

	vt = ippsMalloc_32f(nwx);
	
	ippsFFTInitAlloc_R_32f(&spec, order, IPP_FFT_NODIV_BY_ANY, ippAlgHintNone);
	ippsFFTGetBufSize_R_32f(spec, &bs);
	buf = ippsMalloc_8u(bs);

	ippsFFTFwd_RToPerm_32f_I(other->vecFilterRowR, spec, buf);
	

	printTag(aD,"SliceFirst",af->sliceFirst);
	printTag(aD,"SliceLast",af->sliceLast);
	printTag(aD,"SliceStep",af->sliceStep);
	if(af->sliceStep <0.1){
		printError(aD,"wrong number");
		return false;
	}

	slf = af->sliceFirst;
	sll = af->sliceLast;
	sls = af->sliceStep;

	nbo = af->outputBits;
	fprintf(aD->fp,"Number of bits for the output images: %i.\n",nbo);
	if(nbo != 8 && nbo != 16 && nbo != 32){
		printError(aD, "wrong number");
		return false;
	}

	outfile = af->outputIndexFirst;
	fprintf(aD->fp,"Output index (first): %i.\n",af->outputIndexFirst);
	if(af->outputIndexFirst < 0){
		printError(aD, "wrong number");
		return false;
	}
	fprintf(aD->fp,"Output index (step): %i.\n",af->outputIndexStep);
	if(af->outputIndexStep < 1){
		printError(aD, "wrong number");
		return false;
	}

	
	sprintf(otemplate,"%s/%s%%0%ii%s.%s",af->outputFolder,af->outputPrefix,af->outputNOD,af->outputSuffix,af->outputExtension);

	
	printTagStart(aD,"Slices");

	for(float sl = slf; sl <= sll; sl+=sls){
		printf("Slice %f\n", sl);
		matp = new MatPar();
		fillMatPar(matp,PC,aD,other,sl);
		if(!matp->formCudaMatrixXY())return false;
		if(!readInputVolume(matp, PC, aD, other))return false;

		ny = matp->y_maxi - matp->y_mini+1;

		for(int k=0;k<ny*nz;k++){
			ippsCopy_32f(matp->ivol+k*nx,vt+ext_left,nx);
			ippsMean_32f(matp->ivol+k*nx,ext_num,&sa,ippAlgHintNone);
			if(ext_left>0)ippsSet_32f(sa,vt,ext_left);
			ippsMean_32f(matp->ivol+k*nx+nx-ext_num,ext_num,&sa,ippAlgHintNone);
			if(ext_right>0)ippsSet_32f(sa,vt+ext_left+nx,ext_right);
			ippsFFTFwd_RToPerm_32f_I(vt, spec, buf);
			ippsMulPerm_32f_I(other->vecFilterRowR, vt, nwx);
			ippsFFTInv_PermToR_32f_I(vt, spec, buf);
			ippsCopy_32f(vt+ext_left,matp->ivol+k*nx,nx);
		}
		
		if(!cudaReconstructFDK(matp, PC, aD, other))return false;

		to = new gtiff();
		sprintf(ofile,otemplate,outfile);
		fprintf(aD->fp,"Output file: \"%s\".\n",ofile);
		printf("Output file: \"%s\".\n",ofile);

		switch(nbo){
			case 8:
				outb = to->print_float_8u(ofile,matp->rec,PC->OutputWidth,PC->OutputHeight,af->outputMin,af->outputMax);
				break;
			case 16:
				outb = to->print_float_16u(ofile,matp->rec,PC->OutputWidth,PC->OutputHeight,af->outputMin,af->outputMax);
				break;
			default:
				outb = to->print_float_32f(ofile,matp->rec,PC->OutputWidth,PC->OutputHeight);
				break;
		}

		delete to;

		if(outb){
			printInfo(aD, "The file has been written");
		}else{
			printError(aD, "can not write to the file");
			return false;
		}

		outfile+= af->outputIndexStep;

		delete matp;
	}

	printTagEnd(aD);

	ippsFFTFree_R_32f(spec); spec = NULL;
	ippsFree(buf); buf = NULL;
	ippsFree(vt); vt = NULL;
	
	delete other; other = NULL;
//	delete mp; mp = NULL;

	
	print_global_time_stamp(aD);
	//fclose(aD->fp);
	//aD->fp = NULL;
	return true;

}
