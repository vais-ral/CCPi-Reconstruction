/* ----------------------------------------------------------------------
 * ole2cd.cpp()
 * ----------------------------------------------------------------------
 * 
 * ----------------------------------------------------------------------
 */

#include "defs.h"

void time_to_char(time_stamp *ts, char *c){
	int seci = int(floor(ts->sec));
	int remi = int(floor((ts->sec-float(seci))*1000000.f));
	sprintf(c,"%02i-%02i-%04i %02i:%02i:%02i.%06i", (int)ts->days, (int)ts->months, (int)ts->years, (int)ts->hours, (int)ts->minutes, (int)seci, (int)remi);
}


dir_entry::dir_entry(){
	vtype = 0;
	comm[0]='\0';
	esize=4;
	pnote = true;
	psize = true;
	ptime_mod = true;
	ptime_creation = true;
}

ole2cd::ole2cd(){
	num_child = NULL;
	child_ind = NULL;
	parent =  NULL;
	vvv = NULL;
	emb = NULL;
	emb_ind = NULL;
	MSAT = NULL;
	SSAT = NULL;
	SAT = NULL;
	ff = NULL;
	sectors = NULL;
	cdfID[8] = '\0';
	UID[16] = '\0';
	byte_order = true;
}

ole2cd::~ole2cd(){
	Ipp32u i;
	for(i = 0; i<tot_dir; i++){
		delete p_dir[i];
	}
	if(p_dir != NULL){free(p_dir); p_dir = NULL;}
	if(num_child != NULL){free(num_child); num_child = NULL;}
	if(child_ind != NULL){free(child_ind); child_ind = NULL;}
	if(parent != NULL){free(parent); parent = NULL;}
	if(vvv != NULL){free(vvv); vvv = NULL;}
	if(emb != NULL){free(emb); emb = NULL;}
	if(emb_ind != NULL){free(emb_ind); emb_ind = NULL;}
	if(sectors != NULL){free(sectors); sectors = NULL;}
	if(MSAT != NULL){ippsFree(MSAT); MSAT = NULL;}
	if(SSAT != NULL){ippsFree(SSAT); SSAT = NULL;}
	if(SAT != NULL){ippsFree(SAT); SAT = NULL;}
	if(ff != NULL){fclose(ff); ff = NULL;}
}

bool ole2cd::print_data(allData *aD, char *fraw){
	char fname[MAX_FOLDER], im_name[MAX_FOLDER], fs[MAX_FOLDER];
	int imf, iml, ims, ntot;
	int file_f, file_l, file_s;
	int ipf, v, rnd, nd;
	int rowf, rowl, nx, ny, nt, nc, pss, u, st, nb;
	int nf, st1, st2;
	Ipp32u j;
	
	bool is_ok, bp;
	
	Ipp8u *c;
	Ipp8u *veco;
	gtiff *to;
	xXrm2tif *x2t;
	dir_entry *p;
	
	x2t = aD->hms->x2t;

	printTag(aD,"BitsPerData",x2t->bitsPerData);
		
	if(!(x2t->bitsPerData == 8 || x2t->bitsPerData == 16 || x2t->bitsPerData == 32)){
		printError(aD,"wrong number of bits");
		return false;
	}
	
	imf = __max(0,x2t->imageFirst);
	iml = __min(x2t->imageLast,x2t->NOI-1);

	printTag(aD,"ImageFirst",imf);
	printTag(aD,"ImageLast",iml);
		
	if(iml<imf){
		printError(aD,"the first image > the last image");
		return false;
	}

	ims = __max(1,x2t->imageStep);
	ntot = (iml-imf)/ims+1;
	iml = imf + (ntot-1)*ims;
	ipf = __min(__max(1,x2t->imagesPerFile),ntot);
	file_f = x2t->fileFirst;
	file_s = __max(1,x2t->fileStep);
	file_l = file_f + file_s *((ntot-1)/ipf);


	printTag(aD,"ImageStep",ims);
	printTag(aD,"ImagesTotal",ntot);
	printTag(aD,"ImageFirstUpdated",imf);
	printTag(aD,"ImageLastUpdated",iml);
	printTag(aD,"ImagesPerFile",ipf);	
	printTag(aD,"FileFirst",file_f);
	printTag(aD,"FileLast",file_l);
	printTag(aD,"FileStep",file_s);

	v = file_l;
	rnd = 1;
	while(v>9){
		rnd++;
		v/=10;
	}

	printTag(aD, "NumberOfDigitsForLastFile", rnd);

	nd = x2t->NOD;
	printTag(aD,"NumberOfDigits",nd);
	
	if(nd < 1 || nd >9 || nd<rnd){
		printError(aD,"wrong number");
		return false;
	}
		
	sprintf(fs,"%s/%s%%0%ii%s.%s",fraw,x2t->prefix,nd,x2t->suffix,x2t->extension);
		
	rowf = __max(0,x2t->rowFirst);
	rowl = __min(x2t->height-1,x2t->rowLast);

	printTag(aD,"RowFirst",rowf);
	printTag(aD,"RowLast",rowl);

	if(rowl<rowf){
		printError(aD, "the first row > the last row");
		return false;
	}

	ny = rowl -rowf+1;
	nx = x2t->width;
	nt = nx*ny;

	nb = x2t->bitsPerData/8;
	
	veco = ippsMalloc_8u(nt*ipf*nb);
	ippsZero_8u(veco,nt*ipf*nb);
	if(veco == NULL){
		printError(aD,"can not allocate memory to write the images");
		return false;
	}

	to = new gtiff();
	
	printInfo(aD,"Multiple images to be extracted");

	nf = 0;
	nc = ipf;

	printTagStart(aD,"ImagesAndFiles");

	for(int i = imf; i<=iml; i+=ims){
		if(nf%ipf == 0){
			sprintf(fname,fs,file_f+file_s*(nf/ipf));
			printTag(aD,"File",fname);
		}
		printTag(aD,"Image",i);
		if(i == iml){
			nc = ntot%ipf;
			if(nc == 0) nc = ipf;
		}
		
		sprintf(im_name,"Image%i",i+1);
		is_ok = false;
		
		for(j = 0; j<tot_dir; j++){
			p = this->p_dir[j];
			if(strcmp(p->name,im_name) != 0)continue;
			if(this->parent[j]<1)continue;
			is_ok=true;
			break;
		}
	
		if(!is_ok){
			sprintf(aD->message,"can not find the image %i",i);
			printError(aD);
			return false;
		}

		pss = p->stream_size;
		st1 = (rowf*nb*nx)/size_sector;
		st2 = ((rowf+ny)*nb*nx-1)/size_sector;
		if(i == imf) c =  ippsMalloc_8u(pss+size_sector);
		u = p->secID;
		FSEEK_FUNC(ff,(long long)(u+1)*(long long)(size_sector),SEEK_SET);
		st = 0;

		do{
			if(st >= st1 && st <= st2)fread(c+(size_sector*st),sizeof(char),size_sector,ff);
			st++;	
			pss-=(size_sector);
			if(pss>0) u = SAT[u];
			if(st >= st1 && st <= st2)FSEEK_FUNC(ff,(long long)(u+1)*(long long)(size_sector),SEEK_SET);
		}while(pss>0);

		ippsCopy_8u(c+rowf*nx*nb,veco+(nf%ipf)*nx*ny*nb,nx*ny*nb);

		nf++;
		if(i == iml || nf%ipf == 0){
			switch(x2t->bitsPerData){
				case 8:
					bp = to->print_int_8u(fname,veco,nx,ny*nc);
					break;
				case 16:
					bp = to->print_int_16u(fname,(Ipp16u *)veco,nx,ny*nc);
					break;
				case 32:
					bp = to->print_float_32f(fname,(Ipp32f *)veco,nx,ny*nc);
					break;
				default:
					break;
			}
			if(bp){
				printInfo(aD,"The image file has been created");
			}else{
				printError(aD,"can not write data to the image file");
				return false;
			}
		}
	
	}

	printTagEnd(aD);
	
	delete to;
	ippsFree(veco); veco = NULL;
	ippsFree(c); c = NULL;

	return true;
}



bool ole2cd::print_ref_data(allData *aD){
	char fname[MAX_FOLDER], im_name[MAX_FOLDER];
	xXrm2tif *x2t;
	Ipp32u j;
	int nx, ny, rowf, rowl;
	int pss, u, st;
	Ipp8u *c;
	dir_entry *p;
	bool is_ok, bp;
	gtiff *to;

	x2t = aD->hms->x2t;
	sprintf(fname,"%s/raw/%s.%s",x2t->folder,x2t->flatName,x2t->extension);
	printTag(aD,"ReferenceFile",fname);

	printTag(aD,"BitsPerRef",x2t->bitsPerRef);
		
	if(!(x2t->bitsPerRef == 8 || x2t->bitsPerRef == 16 || x2t->bitsPerRef == 32)){
		printError(aD,"wrong number of bits");
		return false;
	}

	sprintf(im_name,"Image");
	is_ok = false;
	
	for(j = 0; j<tot_dir; j++){
		p = this->p_dir[j];
		if(strcmp(p->name, im_name) != 0)continue;
		if(this->parent[j]<1)continue;
		if(strcmp(p_dir[parent[j]]->name, "ReferenceData") != 0)continue;
		is_ok=true;
		break;
	}
	
	if(!is_ok){
		printError(aD, "can not find the reference image");
		return false;
	}
	
	rowf = __max(0,x2t->rowFirst);
	rowl = __min(x2t->height-1,x2t->rowLast);

	printTag(aD,"RowFirst",rowf);
	printTag(aD,"RowLast",rowl);
	
	if(rowl<rowf){
		printError(aD, "the first row > the last row");
		return false;
	}

	ny = rowl - rowf+1;
	nx = x2t->width;
	
	pss = p->stream_size;
	c =  ippsMalloc_8u(pss+size_sector);
	u = p->secID;
	FSEEK_FUNC(ff,(long long)(u+1)*(long long)(size_sector),SEEK_SET);
	st = 0;
			
	do{
		fread(c+(size_sector*st),sizeof(char),size_sector,ff);
		st++;	
		pss-=(size_sector);
		if(pss>0) u = SAT[u];
		FSEEK_FUNC(ff,(long long)(u+1)*(long long)(size_sector),SEEK_SET);
	}while(pss>0);

	bp = true;
	to = new gtiff();
	switch(x2t->bitsPerRef){
		case 8:
			bp = to->print_int_8u(fname,c+rowf*nx,nx,ny);
			break;
		case 16:
			bp = to->print_int_16u(fname,(Ipp16u *)(c+2*rowf*nx),nx,ny);
			break;
		case 32:
			bp = to->print_float_32f(fname,(Ipp32f *)(c+4*rowf*nx),nx,ny);
			break;
	}
	
	delete to;
	ippsFree(c); c = NULL;

	if(bp){
		printInfo(aD,"The reference file has been created");
	}else{
		printError(aD,"can not write data to the reference file");
		return false;
	}
	return true;
}


void print_sh(FILE *fp, Ipp8u* q, int m, int n, int pe, int vt, char * name){
	//char temp[500];
	//sprintf(temp,"\n<v%s num=\"%%i\">%%f</v%s>",name,name);
	if(vt == 1 || vt == 4){
		fprintf(fp,"%i",*((int *)q));
	}else if(vt == 2){
		fprintf(fp,"%f",*((float *)q));
	}else if(vt == 3){
		fprintf(fp,"%s",q);
	}else if(vt == -1 || vt == -4){
		for(int i=0;i<n;i++)fprintf(fp,"\n<v%s num=\"%i\">%i</v%s>",name,i+m,*((int *)(q+i*pe)),name);
	}else if(vt == -2){
		for(int i=0;i<n;i++)fprintf(fp,"\n<v%s num=\"%i\">%f</v%s>",name,i+m,*((float *)(q+i*pe)),name);
//		for(int i=0;i<n;i++)fprintf(fp,temp,i+m,*((float *)(q+i*pe)));
	}else if(vt == -3){
		for(int i=0;i<n;i++)fprintf(fp,"\n<v%s num=\"%i\">%s</v%s>",name,i+m,(char *)(q+i*pe),name);
	}
}
void print_all(dir_entry *p, ole2cd *od, FILE *fp, Ipp8u *q){
	int uu;
	int k, ssps, u, pss, st;
	int oss, osss, vt, pe, w, m, n, b, r;
	long long lpos;
	
	oss = od->size_sector;
	osss = od->size_short_sector;
	vt = p->vtype;
	pss = p->stream_size;
	u = p->secID;
	pe = p->esize;

	ssps = oss/osss;
	w = pss;
	m = 0;
	n = pss/pe;
	r = 0;
	st = 0;
		
	fprintf(fp,"<%s",p->name);
	if(strlen(p->comm)>0 && p->pnote) fprintf(fp," note=\"%s\"",p->comm);
	if(vt == 0 && p->psize) fprintf(fp," size=\"%i\"",p->stream_size);
	if(p->ptime_creation && p->t_creation.years>1900) fprintf(fp," time_creation=\"%s\"",p->t_creation.cname);
	if(p->ptime_mod && p->t_modification.years>1900) fprintf(fp," time_mod=\"%s\"",p->t_modification.cname);
	fprintf(fp,">");

	if(vt != 0){
		if(pss <= oss){		
			w = pss;
			st = 0;
			for(;;){
				uu = u;
				k = (int)od->start_ss;
				while(uu>=ssps){
					uu -= ssps;
					k = od->SAT[k];
				};
				lpos = (long long)(k+1)*(long long)(oss)+(long long)(osss)*(long long)(uu);
				FSEEK_FUNC(od->ff,lpos,SEEK_SET);
				fread(q+(osss*st),sizeof(Ipp8u),osss,od->ff);
				st++;
				
				w-=osss;
				if(w<=0)break;
				u = od->SSAT[u];
			}

			print_sh(fp, q, m, n, pe, vt, p->name);
		}else{
			FSEEK_FUNC(od->ff,(long long)(u+1)*(long long)(oss),SEEK_SET);
			do{
				fread(q+r, sizeof(Ipp8u), oss, od->ff);
				b = (__min(w,oss)+r)/pe;
				r = __min(w,oss)+ r - b*pe;
				print_sh(fp, q, m, b, pe, vt, p->name);
				if(r>0)	ippsCopy_8u(q+b*pe, q, r);
				m+=b;
				w-=oss;
				if(w>0) u = od->SAT[u];
				FSEEK_FUNC(od->ff,(long long)(u+1)*(long long)(oss),SEEK_SET);
			}while(w>0);
		}
	}
	
	fprintf(fp,"</%s>\n",p->name);
}


void print_top(dir_entry *p, FILE *fp){
	if(p->type == 5){
		fprintf(fp,"<RootEntry");
	}else{
		fprintf(fp,"<%s",p->name);
		if(p->ptime_creation && p->t_creation.years>1900)fprintf(fp," time_creation=\"%s\"",p->t_creation.cname);
		if(p->ptime_mod && p->t_modification.years>1900)fprintf(fp," time_mod=\"%s\"",p->t_modification.cname);
	}
	fprintf(fp,">\n");
}

void print_bottom(dir_entry *p, FILE *fp){
	if(p->type == 5){
		fprintf(fp,"<RootEntry>\n");
	}else{
		fprintf(fp,"<%s>\n",p->name);
	}
}


int ole2cd::pow2(int s){
	if(s<0)return 0;
	return (1<<s);
}

bool ole2cd::read_header(){
	long long lsize_sector;
	Ipp8u header[512];
	Ipp16u s16;
	Ipp32s s, sis;
	Ipp32u k, len_SAT, ss_per_sect, len_SSAT, ids_per_sect;
	Ipp32u sect_head, len_MSAT, num_sect_MSAT;
	Ipp32u dir_per_sect;
	Ipp32u ii, i, j;
	
	fread(header,sizeof(Ipp8u),512,ff);
	
	ippsCopy_8u(header,(Ipp8u *)cdfID,8);
	ippsCopy_8u(header+8,(Ipp8u *)UID,16);
	ippsCopy_8u(header+24,(Ipp8u *)(&revision_number),2);
	ippsCopy_8u(header+26,(Ipp8u *)(&version_number),2);
	ippsCopy_8u(header+28,(Ipp8u *)(&s16),2);//byte order 65534 - little-endian, 65279 - big-endian
	ippsCopy_8u(header+30,(Ipp8u *)(&s16),2);
	size_sector = (Ipp32u)pow2(s16);
	lsize_sector = (long long)size_sector;

	ippsCopy_8u(header+32,(Ipp8u *)(&s16),2);
	size_short_sector = pow2(s16);

	ippsCopy_8u(header+40,(Ipp8u *)(&num_sect_dir),4);
	ippsCopy_8u(header+44,(Ipp8u *)(&num_sec_for_SAT),4);
	ippsCopy_8u(header+48,(Ipp8u *)(&SecID_first_dir),4);
	ippsCopy_8u(header+56,(Ipp8u *)(&min_size_std_stream),4);

	ippsCopy_8u(header+60,(Ipp8u *)(&SecID_first_SSAT),4);
	ippsCopy_8u(header+64,(Ipp8u *)(&num_SSAT),4);
	ippsCopy_8u(header+68,(Ipp8u *)(&SecID_first_MSAT),4);
	ippsCopy_8u(header+72,(Ipp8u *)(&num_MSAT),4);
	
	ids_per_sect = size_sector/4;
	sect_head = (512-76)/4;
	len_MSAT = num_MSAT*ids_per_sect + sect_head;
	num_sect_MSAT = 0;

	MSAT = ippsMalloc_32u((int)len_MSAT);
	
	for(j=0;j<sect_head;j++){
		ippsCopy_8u(header+76+4*j,(Ipp8u*)&s,4);
		if(s<0)break;
		MSAT[j] = (Ipp32u)s;
		num_sect_MSAT++;
	}

	k = sect_head;

	if(num_MSAT>0)FSEEK_FUNC(ff,(long long)(SecID_first_MSAT+1)*lsize_sector,SEEK_SET);

	if(s>=0){
		for(j = 0; j<num_MSAT; j++){
			for(i = 0; i<ids_per_sect; i++){
				fread(&s, sizeof(Ipp32s), 1, ff);
				if(s<0) break;
				MSAT[k] = Ipp32u(s);
				k++;
			}
			if(s<0)break;
			k--;
			if(s>=0)FSEEK_FUNC(ff,(long long)(s+1)*lsize_sector,SEEK_SET);
		}
		num_sect_MSAT = k;
	}

	len_SAT = num_sect_MSAT*ids_per_sect;
	SAT = ippsMalloc_32s((int)len_SAT);

	for(j = 0; j<num_sect_MSAT; j++){
		FSEEK_FUNC(ff,(long long)(MSAT[j]+1)*lsize_sector,SEEK_SET);
		fread(SAT+j*ids_per_sect,sizeof(Ipp32s),ids_per_sect,ff);
	}

	ss_per_sect = size_sector/4;
	
	
	if(num_SSAT>0){
		len_SSAT = num_SSAT*ss_per_sect;
		SSAT = ippsMalloc_32u(len_SSAT);
		sis = (Ipp32s)SecID_first_SSAT;
		for(j=0; j<num_SSAT; j++){
			FSEEK_FUNC(ff,(long long)(sis+1)*lsize_sector,SEEK_SET);
			fread(SSAT+j*ss_per_sect,sizeof(Ipp32u),ss_per_sect,ff);
			sis = SAT[sis];
			if(sis<0)break;
		}
	}

	dir_per_sect = size_sector/128;
	tot_dir = 0;
	s = (Ipp32s)SecID_first_dir;

	for(;;){
		s = SAT[(unsigned int)s];
		tot_dir += dir_per_sect;
		if(s < 0)break;
	}


	p_dir = (dir_entry**)malloc(sizeof(dir_entry*)*tot_dir);
	s = (Ipp32s)SecID_first_dir;
	ii = 0;

	for(;;){
		FSEEK_FUNC(ff,(long long)(s+1)*lsize_sector,SEEK_SET);
		for(j = 0; j < dir_per_sect; j++){
			p_dir[ii] = new dir_entry();			
			read_dir_entry(p_dir[ii]);		
			ii++;
		}
		s = SAT[s];
		if(s<0)break;
	}
	return true;
}

bool ole2cd::sort_dir(){	
	int tk, y1, y2;
	int ti, ew;
	Ipp32s sk;
	bool tq;
	Ipp32u i, m;

	num_child = (int *)malloc(sizeof(int)*tot_dir);
	child_ind = (int *)malloc(sizeof(int)*tot_dir);
	parent = (int *)malloc(sizeof(int)*tot_dir);

	vvv = (int *)malloc(sizeof(int)*tot_dir);
	emb = (int *)malloc(sizeof(int)*tot_dir);
	emb_ind = (int *)malloc(sizeof(int)*tot_dir);

	for(i = 0; i<tot_dir; i++){
		parent[i]=-1;
		num_child[i]=0;
		emb[i]=0;
	}

	for(i = 0; i<tot_dir; i++){
		if(p_dir[i]->type == 0)continue;
		if(p_dir[i]->dirID_root == -1) continue;
		parent[p_dir[i]->dirID_root]=i;		
	}
	
	y1 = 0;
	y2 = 0;
	ew = 0;

	do{
		y2 = y1;
		y1 = 0;
		for(i = 0; i<tot_dir; i++){	
			if(p_dir[i]->type == 0)continue;
			tk = parent[i];
			if(tk == -1) continue;
			sk = p_dir[i]->dirID_left;
			while(sk > -1){
				parent[(int)sk] = tk;
				y1++;
				sk = p_dir[sk]->dirID_left;
			};

			sk = p_dir[i]->dirID_right;
			while(sk>-1){
				parent[(int)sk] = tk;
				y1++;
				sk = p_dir[sk]->dirID_right;
			};
		}
	}while(y1 !=y2);

	for(i = 0; i<tot_dir; i++){
		if(parent[i]>-1)num_child[parent[i]]++;
	}
	
	root_ind=0;
	for(i = 0; i<tot_dir; i++){
		if(p_dir[i]->type == 5)break;
		root_ind++;
	}

	
	for(i = 0; i<tot_dir; i++){
		if(num_child[i] == 0){
			child_ind[i] = -1;
			continue;
		}
		child_ind[i] = ew;
		for(m = 0; m < tot_dir; m++){
			if(parent[m] == (int)i){
				vvv[ew]=m;
				ew++;
			}
		}

		do{
			tq = true;
			for(int k = ew-num_child[i]; k<ew-1; k++){
				if(strncmp(p_dir[vvv[k]]->name,p_dir[vvv[k+1]]->name,32)>0){
					ti = vvv[k];
					vvv[k]=vvv[k+1];
					vvv[k+1] = ti;
					tq = false;
				}
			}
		}while(!tq);
	}

	this->start_ss = p_dir[root_ind]->secID;
	return true;
}





bool ole2cd::fill_MainParam(allData *aD){	
	xXrm2tif *xt; 
	this->start_ss = p_dir[root_ind]->secID;
	long long oss, osss, offset;
	Ipp32s xu, xk, ssps, pss;
	Ipp32u i;
	int v1, tag;

	xt = aD->hms->x2t;
	
	oss = (long long)size_sector;
	osss = (long long)size_short_sector;
	ssps = (Ipp32s) (oss/osss);

	for(i = 0; i<tot_dir; i++){
		tag = XRM_UNKNOWN;
		if(strcmp("NoOfImages",p_dir[i]->name) == 0) tag = XRM_NOI;
		if(strcmp("ImageHeight",p_dir[i]->name) == 0) tag = XRM_HEIGHT;
		if(strcmp("ImageWidth",p_dir[i]->name) == 0) tag = XRM_WIDTH;
		if(strcmp("PixelSize",p_dir[i]->name) == 0) tag = XRM_PIXEL_SIZE;
		if(strcmp("DataType",p_dir[i]->name) == 0) tag = XRM_DATA_TYPE;
		if(strcmp("DtoRADistance",p_dir[i]->name) == 0) tag = XRM_D2A;
		if(strcmp("StoRADistance",p_dir[i]->name) == 0) tag = XRM_S2A;
		if(tag == XRM_UNKNOWN) continue;

		xu = p_dir[i]->secID;
		xk = start_ss;
		if(tag == XRM_D2A || tag == XRM_S2A){
			pss = p_dir[i]->stream_size;

			if(pss<= (Ipp32s)size_sector){
				while(xu >= ssps){
					xu -= ssps;
					xk = this->SAT[xk];
				};
				offset = (long long)(xk+1)*oss + osss*(long long)(xu);
			}else{
				offset = (long long)(xu+1)*oss;
			}
			FSEEK_FUNC(this->ff, offset, SEEK_SET);
			if(parent[i]>-1){
				if(strcmp("ImageInfo",p_dir[parent[i]]->name) == 0){
					if(tag == XRM_D2A){
						read_float(&(xt->d2a));
					}else{
						read_float(&(xt->s2a));
					}
				}
			}	
		}else{
			while(xu>=ssps){
				xu -= ssps;
				xk = this->SAT[xk];
			};
			offset = (long long)(xk+1)*oss+osss*(long long)(xu);
			FSEEK_FUNC(this->ff,offset,SEEK_SET);
			switch(tag){
				case XRM_NOI:
					if(parent[i]>-1){
						if(strcmp("ImageInfo",p_dir[parent[i]]->name) == 0)	read_long(&(xt->NOI));
					}
					break;
				case XRM_HEIGHT:
					read_long(&(xt->height));
					break;
				case XRM_WIDTH:
					read_long(&(xt->width));
					break;
				case XRM_PIXEL_SIZE:
					if(parent[i]>-1){
						if(strcmp("ImageInfo",p_dir[parent[i]]->name) == 0)	read_float(&(xt->pixelSize));
					}
					break;
				case XRM_DATA_TYPE:
					if(parent[i]>-1){
						if(strcmp("ImageInfo",p_dir[parent[i]]->name) == 0)	read_long(&(xt->bitsPerData));
						if(strcmp("ReferenceData",p_dir[parent[i]]->name) == 0) read_long(&(xt->bitsPerRef));
					}	
					break;
			}

		}
	}

	v1 = xt->bitsPerData;
	
	if(v1 == 3){
		v1 = 8;
	}else if(v1 == 5){
		v1 = 16;
	}else if(v1 == 10){
		v1 = 32;
	}
	xt->bitsPerData = v1;

	v1 = xt->bitsPerRef;
	
	if(v1 == 3){
		v1 = 8;
	}else if(v1 == 5){
		v1 = 16;
	}else if(v1 == 10){
		v1 = 32;
	}
	xt->bitsPerRef = v1;
	
	return true;
}

bool ole2cd::print_MainParam(allData *aD, char *finfo){
	char filen[MAX_FOLDER];
	FILE *fp;
	xXrm2tif *xt;

	xt = aD->hms->x2t;

	sprintf(filen,"%s/main_param.xml",finfo);
	printTag(aD,"File",filen);
	
	fp = fopen(filen,"w");
	if(fp == NULL){
		printError(aD,"can not open the file");
		return false;
	}

	fprintf(fp, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
	fprintf(fp,"<HMxml>\n<MainParam>\n");
	fprintf(fp,"<ImageWidth>%i</ImageWidth>\n",xt->width);
	fprintf(fp,"<ImageHeight>%i</ImageHeight>\n",xt->height);
	fprintf(fp,"<NumberOfImages>%i</NumberOfImages>\n",xt->NOI);
	fprintf(fp,"<PixelSize>%f</PixelSize>\n",xt->pixelSize);
	fprintf(fp,"<DistanceSourceToAxis>%f</DistanceSourceToAxis>\n",xt->s2a);
	fprintf(fp,"<DistanceDetectorToAxis>%f</DistanceDetectorToAxis>\n",xt->d2a);
	fprintf(fp,"<BitsPerData>%i</BitsPerData>\n",xt->bitsPerData);
	if(xt->bitsPerRef >0)fprintf(fp,"<BitsPerRef>%i</BitsPerRef>\n",xt->bitsPerRef);
	fprintf(fp,"</MainParam>\n</HMxml>\n");
	fclose(fp); fp = NULL;

	return true;
}


bool ole2cd::read_types(allData *aD){
	char *buff, *buff2, *b1, *b2, *b3, *ctc, *pdest;
	char cxml[5] = "?xml";
	int p1, p2, w1, w2, w3;
	bool tb, tc;
	Ipp32u i;
		
	FILE *ft;

	tb = true;

	printTag(aD,"TagFile",aD->hms->x2t->tagFile);
		
	ft = fopen(aD->hms->x2t->tagFile,"r");
	if(ft == NULL){
		printError(aD,"can not open the file");
		return false;
	}
	printInfo(aD,"The file has been opened");

	buff = (char *)malloc(line_width*sizeof(char));
	buff2 = (char *)malloc(line_width*sizeof(char));

	b1 = (char *)malloc(line_width*sizeof(char));
	b2 = (char *)malloc(line_width*sizeof(char));
	b3 = (char *)malloc(line_width*sizeof(char));
	ctc = (char *)malloc(line_width*sizeof(char));

	tc = true;

	while(fgets(buff,line_width,ft) != NULL){
		if(int(strlen(buff))<2)continue;
		pdest = strchr(buff, '<');
		if(pdest == NULL)continue;
		p1 = (int)(pdest - buff);
		pdest = strrchr(buff, '>');
		if(pdest == NULL)continue;
		p2 = (int)(pdest - buff);
		if(p2-p1<3)continue;
		if(strncpy(buff2, buff+p1+1, p2-p1-1) == NULL)continue;
		buff2[p2-p1-1] = '\0';
		pdest = strrchr(buff2, '>');
		b2[0] = '\0';
		b3[0] = '\0';
		w1=0;
		w2=0;
		w3=0;
		if(pdest == NULL){
			w1=p2-1;
			w2=0;
			w3=0;
			strncpy(b1, buff2, w1+1);
			if(strncmp(b1, cxml ,4) == 0) continue;
			b1[w1] = '\0';
			b2[0] = '\0';
			b3[0] = '\0';
			if(tc){
				tc = false;
				if(strcmp(b1,"int") == 0){
					strcpy(ctc,"/int");
				}else if(strcmp(b1,"float") == 0){
					strcpy(ctc,"/float");
				}else if(strcmp(b1,"char") == 0){
					strcpy(ctc,"/char");
				}else if(strcmp(b1,"bool") == 0){
					strcpy(ctc,"/bool");
				}else if(strcmp(b1,"array") == 0){
					strcpy(ctc,"/array");
				}else if(strcmp(b1,"ignore") == 0){
					strcpy(ctc,"/ignore");
				}else if(strcmp(b1,"exclude") == 0){
					strcpy(ctc,"/exclude");
				}else{
					tc = true;
				}
			}else{
				if(strcmp(ctc,"/int") == 0 || strcmp(ctc,"/float") == 0 || strcmp(ctc,"/char") == 0 || strcmp(ctc,"/bool") == 0)tc = true;
				if(strcmp(ctc,"/array") == 0 || strcmp(ctc,"/ignore") == 0 || strcmp(ctc,"/exclude") == 0)tc = true;
			}
		}else if(!tc){
			p2 = (int)(pdest - buff2);
			w1 = p2;
			strncpy(b1, buff2, w1);
			b1[w1] = '\0';
			p2 = int(strlen(buff2))-w1;
			strncpy(buff, buff2+w1+1, p2-1);
			buff[p2-1] = '\0';
			pdest = strrchr(buff, '/');
			p2 = (int)(pdest - buff);
			w3 = int(strlen(buff))-p2-1;
			strncpy(b3, buff+p2+1, w3);
			b3[w3] = '\0';
			w2 = p2-1;
			strncpy(b2, buff, w2);
			b2[w2] = '\0';
			if(strcmp(b1,b3) != 0)continue;
			tb = false;
			if(strcmp(b1, "note") == 0){
				for(i = 0; i<tot_dir; i++){
					p_dir[i]->pnote = false;
				}
			}
			if(strcmp(b1,"size") == 0){
				for(i = 0; i<tot_dir; i++){
					p_dir[i]->psize = false;
				}
			}
			if(strcmp(b1,"time_mod") == 0){
				for(i = 0; i<tot_dir; i++){
					p_dir[i]->ptime_mod = false;
				}
			}
			if(strcmp(b1,"time_creation") == 0){
				for(i = 0; i<tot_dir; i++){
					p_dir[i]->ptime_creation = false;
				}
			}

			for(i = 0; i<tot_dir; i++){
				if(strcmp(b1,p_dir[i]->name) == 0){	
					if(strcmp(ctc,"/int") == 0){
						p_dir[i]->vtype = 1;
						tb = true;
					}else if(strcmp(ctc,"/float") == 0){
						p_dir[i]->vtype = 2;
						tb = true;
					}else if(strcmp(ctc,"/char") == 0){
						p_dir[i]->vtype = 3;
						tb = true;
					}else if(strcmp(ctc,"/bool") == 0){
						p_dir[i]->vtype = 4;
						tb = true;
					}else if(strcmp(ctc,"/ignore") == 0){
				//		p_dir[i]->vtype = 5;
						tb = true;
					}else if(strcmp(ctc,"/array") == 0){
						if(p_dir[i]->vtype>0) p_dir[i]->vtype *= -1;
						if(atoi(b2) !=0) p_dir[i]->esize = atoi(b2);
						tb = true;
					}
					if(tb){
						if(w2>0 && strcmp(ctc,"/array") != 0)strcpy(p_dir[i]->comm,b2);
					}
				}
			}

			pdest = strrchr(b1, '#');
			if(pdest != NULL){
				p2 = (int)(pdest - b1);
				strncpy(b2,b1,p2);
				b2[p2]='\0';

				for(i = 0; i<tot_dir; i++){
					if(strncmp(b2,p_dir[i]->name,p2) != 0)continue;
					w2 = int(strlen(p_dir[i]->name))-p2;
					strncpy(b1,p_dir[i]->name+p2,w2);
					b1[w2]='\0';
					if(atoi(b1) == 0)continue;
					p_dir[i]->vtype = 5;
				}
			}	
		}
	}
	free(buff); buff = NULL;
	free(buff2); buff2 = NULL;
	free(b1); b1 = NULL;
	free(b2); b2 = NULL;
	free(b3); b3 = NULL;
	free(ctc); ctc= NULL;
	fclose(ft); ft = NULL;
	printInfo(aD,"The tag file has been closed");
	return true;
}

bool ole2cd::print_xml(allData *aD, char *sw){
	char cfull[MAX_FOLDER], cshort[MAX_FOLDER], t1[MAX_FOLDER], t2[MAX_FOLDER];
	int tr, last, iemb;
	int level, ind, m, mw;
	int i_s[10], i_l[10];
	FILE *fu, *fs, *ft;
	Ipp32u i;
	bool tb;
	dir_entry *p;
	Ipp8u *q;
	

	tr = root_ind;
	last = tr;
	iemb = 0;

	sprintf(cfull,"%s/info_full.xml",sw);
	sprintf(cshort,"%s/info_short.xml",sw);

	printTag(aD,"InfoFullFile",cfull);
	fu = fopen(cfull,"w");
	if(fu == NULL){
		printError(aD,"can not open the file");
		return false;
	}

	printTag(aD,"InfoShortFile",cshort);
	fs = fopen(cshort,"w");
	if(fs == NULL){
		printError(aD,"can not open the file");
		return false;
	}

	tb = true;

	q = ippsMalloc_8u(2*this->size_sector);

	for(;;){
		if(tb){
			emb_ind[iemb] = tr;		
			if(p_dir[tr]->type == 2){
				if(p_dir[tr]->vtype <0){
					print_all(p_dir[tr],this,fu,q);
				}else if(p_dir[tr]->vtype <5){
					print_all(p_dir[tr],this,fs,q);
					print_all(p_dir[tr],this,fu,q);
				}
					
				if(emb[iemb] == num_child[last]-1){
					emb[iemb]=0;
					tr=last;
					last=emb_ind[iemb-2];
					tb=false;
					iemb--;
					continue;
				}
				emb[iemb]++;
				tr = vvv[child_ind[last]+emb[iemb]];
			
			}else{
				if(p_dir[tr]->vtype != 5){
					print_top(p_dir[tr],fu);
					print_top(p_dir[tr],fs);
				}
				if(num_child[tr] >0){
					if(p_dir[tr]->vtype == 5){
						if(iemb == 0)break;
						tb = true;
						emb[iemb]++;
						if(emb[iemb] == num_child[last]){
							tb=false;
							emb_ind[iemb] = 0;
							tr = last;
							iemb--;
							last = emb_ind[iemb];
							continue;
						}
						tr = vvv[child_ind[last]+emb[iemb]];
					}else{
						iemb++;
						last = tr;
						tr = vvv[child_ind[tr]];
						tb = true;
					}
					continue;
				}
				if(p_dir[tr]->vtype != 5){
					print_bottom(p_dir[tr],fu);
					print_bottom(p_dir[tr],fs);
				}
				emb[iemb]++;
				tb = false;
			}
			continue;
		}else{
			print_bottom(p_dir[tr],fu);
			print_bottom(p_dir[tr],fs);
			if(iemb == 0)break;
			emb[iemb]++;
			tr=last;
			last=emb_ind[iemb-1];
			if(emb[iemb] == num_child[last]){
				tb=false;
				iemb--;
				continue;
			}	
			tr = vvv[child_ind[last]+emb[iemb]];
			tb=true;
		}
	}

	fclose(fu); fu = NULL;
	fclose(fs); fs = NULL;

	
	for(i = 0; i<tot_dir; i++){
		p = this->p_dir[i];
		if(p->vtype>-1 || p->vtype <-4)continue;
		level = 0;
		ind = i;
		t2[0] = '\0';
		m = 0;
		mw = 0;
		do{
			m = int(strlen(p_dir[ind]->name));
			strncpy(t2+mw,p_dir[ind]->name,m);
			i_s[level] = mw;
			i_l[level] = m+1;
			mw +=(m+1);
			t2[mw-1] = '\0';
			level++;			
			ind = parent[ind];
		}while(ind>-1);
		m = sprintf(t1,"%s/",sw);
		for(int j = level-2;j>=0;j--){
			strncpy(t1+m,t2+i_s[j],i_l[j]-1);
			m+=(i_l[j]);
			if(j>0)t1[m-1] = '_';
		}
		t1[m-1]='.';
		t1[m]='x';
		t1[m+1]='m';
		t1[m+2]='l';
		t1[m+3]='\0';
				

		printTag(aD,"File",t1);
		ft = fopen(t1,"w");
		if(ft == NULL){
			printError(aD,"can not write to the file");
			return false;
		}

		fprintf(ft,"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
		for(int j = level-1; j>=1; j--){
			fprintf(ft,"<%s>\n", t2+i_s[j]);
		}

		print_all(p,this,ft,q);

		for(int j = 1; j<=level-1;j++){
			fprintf(ft,"</%s>\n",t2+i_s[j]);
		}
		
		fclose(ft); ft = NULL;
		
	}

	ippsFree(q); q = NULL;

	return true;
}



void ole2cd::read_dir_entry(dir_entry *d){
	char p[1];
	fread(q,sizeof(Ipp8u),128,ff);
	int u=0;
	for(int k=0; k<32; k++){
		wctomb(p,*(wchar_t *)(q+2*k));
		if(p[0] == ' '){
			u++;
			continue;
		}
		d->name[k-u]=p[0];
	}
	
	ippsCopy_8u(q+64,(Ipp8u *)(&(d->size_area)),2);
	d->type = q[66];
	d->node_colour = q[67];

	ippsCopy_8u(q+68,(Ipp8u *)(&(d->dirID_left)),4);
	ippsCopy_8u(q+72,(Ipp8u *)(&(d->dirID_right)),4);
	ippsCopy_8u(q+76,(Ipp8u *)(&(d->dirID_root)),4);
	
	d->t_creation.read_time(q+100);
	d->t_modification.read_time(q+108);
	ippsCopy_8u(q+116,(Ipp8u *)(&(d->secID)),4);
	ippsCopy_8u(q+120,(Ipp8u *)(&(d->stream_size)),4);
}

void dir_entry::print_dir_entry(FILE *fp){
	
	//fwprintf(fp,L"Entry name: %ls\n",name);
	fprintf(fp,"Entry name: %s\n",name);
	fprintf(fp,"Area size: %i\n",size_area);
	fprintf(fp,"Type: %i ",type);
	switch(type){
		case 0:
			fprintf(fp,"(Empty)\n");
			break;
		case 1:
			fprintf(fp,"(User storage)\n");
			break;
		case 2:
			fprintf(fp,"(User stream)\n");
			break;
		case 3:
			fprintf(fp,"(LockBytes (unknown))\n");
			break;
		case 4:
			fprintf(fp,"(Property (unknown))\n");
			break;
		case 5:
			fprintf(fp,"(Root storage)\n");
			break;
		default:
			fprintf(fp,"(Unknown)\n");
			break;
	};

	fprintf(fp,"Node colour: %i ",node_colour);
	switch(node_colour){
		case 0:
			fprintf(fp,"(red)\n");
			break;
		case 1:
			fprintf(fp,"(black)\n");
			break;
		default:
			fprintf(fp,"(unknown)\n");
			break;
	};

	fprintf(fp,"DirID of the left node: %i",dirID_left);
	if(this->dirID_left == -1){
		fprintf(fp," (no)\n");
	}else{
		fprintf(fp,"\n");
	}
	fprintf(fp,"DirID of the right node: %i",dirID_right);
	if(this->dirID_right == -1){
		fprintf(fp," (no)\n");
	}else{
		fprintf(fp,"\n");
	}
	fprintf(fp,"DirID of the root node: %i",dirID_root);
	if(this->dirID_root == -1){
		fprintf(fp," (no)\n");
	}else{
		fprintf(fp,"\n");
	}
	fprintf(fp,"Time of creation:%4i-%2i-%2i %2i:%2i:%f\n", (int)t_creation.years, (int)t_creation.months, (int)t_creation.days, (int)t_creation.hours, (int)t_creation.minutes, t_creation.sec);
	fprintf(fp,"Time of modification:%4i-%2i-%2i %2i:%2i:%f\n", (int)t_modification.years, (int)t_modification.months, (int)t_modification.days, (int)t_modification.hours, (int)t_modification.minutes, t_modification.sec);
	
	
	fprintf(fp,"SecID: %i\n",secID);
	fprintf(fp,"Total size: %i (bytes)\n\n",stream_size);
	
}


void time_stamp::read_time(Ipp8u *q){
	unsigned long long t, tt;
	unsigned int m[12] = {31,28,31,30,31,30,31,31,30,31,30,31};
	
	ippsCopy_8u(q,(Ipp8u *)&t,8);
	tt = t/600000000;
	sec = float(double(t-600000000*tt)/600000000.0);
	minutes = (int)(tt%60);
	t = tt/60;
	hours = (int)(t%24);
	tt = t/24;
	years = 1601;
	unsigned int s;
	for(;;){
		s = 365;
		if(years%4 == 0 && (years%100 != 0 || years%400 == 0))s=366;
		if(tt < s){
			break;
		}
		tt -= s;
		years++;
	}
	if(s == 366) m[1]++;
	months = 1;
	for(int i=0; i<12; i++){
		if(tt<m[i])break;
		months++;
		tt-=m[i];
	}
	days = (int)tt+1;

	int seci = int(floor(sec));
	int remi = int(floor((sec-float(seci))*1000000.f));
	sprintf(cname,"%02i-%02i-%04i %02i:%02i:%02i.%06i",(int)days, (int)months, (int)years, (int)hours, (int)minutes, (int)seci, (int)remi);
}

bool ole2cd::read_short(int *s){
	int a[2];
	a[0] = getc(ff);
	a[1] = getc(ff);

	if(byte_order){
		*s = a[0]+256*a[1];
	}else{
		*s = a[1]+256*a[0];
	}
	return true;
}


bool ole2cd::read_long(int *s){
	fread(s,sizeof(int),1,ff);
	return true;
}


bool ole2cd::read_bool(int *s){
	int ss;
	fread(&ss,sizeof(int),1,ff);
	if(ss == 0){
		*s=0;
	}else{
		*s=1;
	}
	return true;
}


bool ole2cd::read_float(float *f){
	fread(f,sizeof(float),1,ff);
	return true;
}

