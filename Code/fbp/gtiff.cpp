/* ----------------------------------------------------------------------
 * gtiff.cpp()
 * ----------------------------------------------------------------------
 * Read and write a small set of tiff file formats.
 *
 * TODO: Merge xtiff and gtiff classes. Interface to libtiff.
 * Valeriy: gtiff will be removed later. Now it is used for "cone" functions only.
 * ----------------------------------------------------------------------
 */

#include "defs.h"


gtiff::gtiff(){
	file_size = 0;
	sprintf(err_message,"no errors");
	vect = NULL;
	ff = NULL;
	StripOffsets = NULL;
	StripByteCounts = NULL;
	isII = true;
	ImageWidth=0;
	ImageLength=0;
	BitsPerSample=16;
	SamplesPerPixel = 1;
	Compression=1;
	PhotometricInterpretation=1;
	OffsetStripOffsets=0;
	RowsPerStrip=1;
	OffsetStripByteCounts=0;
	OffsetXResolution=0;
	OffsetYResolution=0;
	XResolution[0]=1;
	XResolution[1]=1;
	YResolution[0]=1;
	YResolution[1]=1;
	ResolutionUnit=1;
	SampleFormat=1;
	num_entries = 13;
}

bool gtiff::find_next_ifd(Ipp32u *off){
	Ipp16u num_of_dir;
	read_short(&num_of_dir);
	num_entries+= Ipp32u (num_of_dir);

	for(Ipp16u i=0;i<num_of_dir;i++)
		read_entry();

	read_long(off);
	return true;
}


bool gtiff::ffind_next_ifd(Ipp32u *off){
	Ipp16u num_of_dir;
	fread_short(&num_of_dir);
	num_entries += Ipp32u (num_of_dir);

	for(Ipp16u i=0;i<num_of_dir;i++)
		fread_entry();

	fread_long(off);
	return true;
}


void gtiff::setError(char *err_string){
	sprintf(err_message,err_string);
	printf("Error: %s.\n",err_message);
}

bool gtiff::getInfo(char *fileName, unsigned int *nx, unsigned int *ny, unsigned int *nb){
	if(ff != NULL){
		setError("another file is open");
		return false;
	}
	ff = fopen(fileName,"rb");
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	long new_file_size;
	fseek(ff, 0, SEEK_END);
	new_file_size = ftell(ff);
	
	
	if(new_file_size <=10){	
		sprintf(err_message,"file size (%li)",new_file_size);
		return false;
	}

	if(vect != NULL || new_file_size != file_size){
		ippsFree(vect); vect = NULL;
		file_size = new_file_size;
	}

	vect = ippsMalloc_8u(file_size);
	fseek(ff,0,SEEK_SET);
	fread(vect,sizeof(Ipp8u),file_size,ff);
	if(!fillTags()){
		ippsFree(vect); vect = NULL;
		return false;
	}
	
	*nx = ImageWidth;
	*ny = ImageLength;
	*nb = BitsPerSample;

	return true;

}

bool gtiff::fread_short(Ipp16u *a){
	fread(a,sizeof(Ipp16u),1,ff);
	//*a = *((Ipp16u*)(vect+cur_pos));
	if(isII == false)ippsSwapBytes_16u_I(a,1);
	cur_pos+=2;
	return true;
}

bool gtiff::read_short(Ipp16u *a){
	*a = *((Ipp16u*)(vect+cur_pos));
	if(isII == false)ippsSwapBytes_16u_I(a,1);
	cur_pos+=2;
	return true;
}

bool gtiff::fread_long(Ipp32u *a){
	fread(a,sizeof(Ipp32u),1,ff);
//	*a = *((Ipp32u*)(vect+cur_pos));
	if(isII == false)ippsSwapBytes_32u_I(a,1);
	cur_pos+=4;
	return true;
}


bool gtiff::read_long(Ipp32u *a){
	*a = *((Ipp32u*)(vect+cur_pos));
	if(isII == false)ippsSwapBytes_32u_I(a,1);
	cur_pos+=4;
	return true;
}


bool gtiff::fillTags(){
	
	if(vect[0] != vect[1]){
		sprintf(ctemp,"this is not TIFF, first two byte are %i, %i",vect[0],vect[1]);
		setError(ctemp);
		return false;
	}
	if(vect[0] == 73){
		isII = true;
	}else  if(vect[1] == 77){
		isII = false;
	}else{
		sprintf(ctemp,"this is not TIFF, first two byte are %i, %i",vect[0],vect[1]);
		setError(ctemp);
		return false;
	}
	cur_pos = 2;

	unsigned short sh;
	read_short(&sh);
	
	if(sh != 42){
		sprintf(ctemp,"this is not TIFF, 3rd and 4th bytes give %i",sh);
		setError(ctemp);
	}
	Ipp32u off;
	read_long(&off);
	cur_pos = off;

	do{
		find_next_ifd(&off);
	}while(off >0);
	
	StripByteCountsPerSample=BitsPerSample/8;

	if(StripOffsets != NULL){
		ippsFree(StripOffsets); StripOffsets = NULL;
	}
	if(StripByteCounts != NULL){
		ippsFree(StripByteCounts); StripByteCounts = NULL;
	}
	StripOffsets = ippsMalloc_32u(StripsPerImage);
	StripByteCounts = ippsMalloc_32u(StripsPerImage);

	if(StripsPerImage == 1){
		StripOffsets[0] = OffsetStripOffsets;
		StripByteCounts[0] = StripByteCountsPerSample*ImageWidth*ImageLength;
	}else{
		cur_pos = OffsetStripOffsets;
		for(Ipp32u i=0;i<StripsPerImage;i++)
			read_long(StripOffsets + i);
	
		cur_pos = OffsetStripByteCounts;
		for(Ipp32u i=0;i<StripsPerImage;i++){
			read_long(StripByteCounts +i);
		}
	}
	
	return true;
}

bool gtiff::ffillTags(){
	Ipp8u v0, v1;
	v0 = (Ipp8u)fgetc(ff);
	v1 = (Ipp8u)fgetc(ff);
	if(v0 != v1){
		sprintf(ctemp,"this is not TIFF, first two byte are %i, %i",v0,v1);
		setError(ctemp);
		return false;
	}
	if(v0 == 73){
		isII = true;
	}else  if(v1 == 77){
		isII = false;
	}else{
		sprintf(ctemp,"this is not TIFF, first two byte are %i, %i",v0,v1);
		setError(ctemp);
		return false;
	}
	cur_pos = 2;

	unsigned short sh;
	fread_short(&sh);
	
	if(sh != 42){
		sprintf(ctemp,"this is not TIFF, 3rd and 4th bytes give %i",sh);
		setError(ctemp);
	}
	Ipp32u off;
	fread_long(&off);
	cur_pos = off;

	do{
		ffind_next_ifd(&off);
	}while(off >0);
	
	StripByteCountsPerSample=BitsPerSample/8;

	StripOffsets = ippsMalloc_32u(StripsPerImage);
	StripByteCounts = ippsMalloc_32u(StripsPerImage);

	if(StripsPerImage == 1){
		StripOffsets[0] = OffsetStripOffsets;
		StripByteCounts[0] = StripByteCountsPerSample*ImageWidth*ImageLength;
	}else{
		cur_pos = OffsetStripOffsets;
		for(Ipp32u i=0;i<StripsPerImage;i++)
			fread_long(StripOffsets + i);
	
		cur_pos = OffsetStripByteCounts;
		for(Ipp32u i=0;i<StripsPerImage;i++){
			fread_long(StripByteCounts +i);
		}
	}
	
	return true;
}




bool gtiff::read_entry(){
	Ipp16u tag, type;
	Ipp32u count, value;
	
	read_short(&tag);
	read_short(&type);
	read_long(&count);
	read_long(&value);

//	printf("Tag: %i, type: %i, count: %i, value: %i\n",tag,type,count,value);

	switch(tag){
		case 256:
			ImageWidth=value;
			if(ImageWidth >65536) ImageWidth/=65536;
			break;
		case 257:
			ImageLength =value;
			if(ImageLength >65536) ImageLength/=65536;
			break;
		case 258:
			if(value > 256){
				BitsPerSample=value/65536;
			}else{
				BitsPerSample=value;
			}
			break;
		case 259:
			Compression=value;
			break;
		case 262:
			PhotometricInterpretation=value;
			break;
		case 273:
			StripsPerImage=count;
			OffsetStripOffsets=value;
			break;
		case 278:
			RowsPerStrip=value;
			break;
		case 279:
			OffsetStripByteCounts=value;
			break;
		case 282:
			OffsetXResolution=value;
			break;
		case 283:
			OffsetYResolution=value;
			break;
		case 296:
			ResolutionUnit=value;
			break;
		default:
		break;
	}
	return true;
	
}


bool gtiff::fread_entry(){
	Ipp16u tag, type;
	Ipp32u count, value;
	
	fread_short(&tag);
	fread_short(&type);
	fread_long(&count);
	fread_long(&value);

//	printf("Tag: %i, type: %i, count: %i, value: %i\n",tag,type,count,value);

	switch(tag){
		case 256:
			ImageWidth=value;
			if(ImageWidth >65536) ImageWidth/=65536;
			break;
		case 257:
			ImageLength =value;
			if(ImageLength >65536) ImageLength/=65536;
			break;
		case 258:
			if(value > 256){
				BitsPerSample=value/65536;
			}else{
				BitsPerSample=value;
			}
			break;
		case 259:
			Compression=value;
			break;
		case 262:
			PhotometricInterpretation=value;
			break;
		case 273:
			StripsPerImage=count;
			OffsetStripOffsets=value;
			break;
		case 278:
			RowsPerStrip=value;
			break;
		case 279:
			OffsetStripByteCounts=value;
			break;
		case 282:
			OffsetXResolution=value;
			break;
		case 283:
			OffsetYResolution=value;
			break;
		case 296:
			ResolutionUnit=value;
			break;
		default:
		break;
	}
	return true;
	
}


gtiff::~gtiff(){
	if(ff != NULL) fclose(ff);
	if(StripOffsets != NULL){
		ippsFree(StripOffsets); StripOffsets = NULL;
	}
	if(StripByteCounts != NULL){
		ippsFree(StripByteCounts); StripByteCounts = NULL;
	}
	if(vect != NULL){
		ippsFree(vect); vect = NULL;
	}
}


bool gtiff::read_float(char *fileName, Ipp32f *vec){
	int nt, buff_size, iup, rem_buff_size;

	if(ff != NULL)fclose(ff);
	ff = fopen(fileName,"rb");
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	long new_file_size;
	fseek(ff, 0, SEEK_END);
	new_file_size = ftell(ff);
	
	if(new_file_size <=10){	
		sprintf(err_message,"file size (%li)",new_file_size);
		return false;
	}

	if(vect != NULL || new_file_size != file_size){
		ippsFree(vect); vect = NULL;
		file_size = new_file_size;
	}

	vect = ippsMalloc_8u(file_size);
	fseek(ff,0,SEEK_SET);
	fread(vect,sizeof(Ipp8u),file_size,ff);
	if(!fillTags()){
		ippsFree(vect); vect = NULL;
		return false;
	}
	
	nt=ImageWidth*ImageLength;

	buff_size=8*StripByteCounts[0]/BitsPerSample;
	iup = (int)StripsPerImage;
	if(iup == 1) buff_size = 8*this->OffsetStripByteCounts/BitsPerSample;
	rem_buff_size = nt%buff_size;
	

	if(rem_buff_size == 0){
		rem_buff_size = buff_size;
		if(nt>buff_size)iup--;
	}

	switch(BitsPerSample){
		case 8:
			for(int i=0;i<iup-1;i++){
				ippsConvert_8u32f(vect+StripOffsets[i], vec+i*buff_size, buff_size);
			}
			ippsConvert_8u32f(vect+StripOffsets[iup-1], vec+(iup-1)*buff_size, rem_buff_size);
			break;
		case 16:
			for(int i=0;i<iup-1;i++){
				if(!isII)ippsSwapBytes_16u_I((Ipp16u*)(vect+StripOffsets[i]), buff_size);
				ippsConvert_16u32f((Ipp16u*)(vect+StripOffsets[i]), vec+i*buff_size, buff_size);
			}
			if(!isII)ippsSwapBytes_16u_I((Ipp16u*)(vect+StripOffsets[iup-1]), buff_size);
			ippsConvert_16u32f((Ipp16u*)(vect+StripOffsets[iup-1]), vec+(iup-1)*buff_size, rem_buff_size);
			break;
		case 32:
			for(int i=0;i<iup-1;i++){
				if(!isII)ippsSwapBytes_32u_I((Ipp32u*)(vect+StripOffsets[i]), buff_size);
				ippsCopy_32f((Ipp32f*)(vect+StripOffsets[i]), vec+i*buff_size, buff_size);
			}
			if(!isII)ippsSwapBytes_32u_I((Ipp32u*)(vect+StripOffsets[iup-1]), buff_size);
			ippsCopy_32f((Ipp32f*)(vect+StripOffsets[iup-1]), vec+(iup-1)*buff_size, rem_buff_size);
			break;

		default:
			ippsZero_32f(vec,nt);
			break;
	}
	
	return true;
}


bool gtiff::read_chunk(char *fileName, Ipp8u *vec, int rmin, int rmax){
	if(ff != NULL)fclose(ff);
	ff = fopen(fileName,"rb");
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	
	if(!ffillTags()){
		return false;
	}

	int rs, s1, s2, r1, r2, nbytes;
	nbytes = BitsPerSample/8;
	rs = StripByteCounts[0]/(nbytes*ImageWidth);

	r1 = rmin%rs;
	s1 = rmin/rs;
	r2 = rmax%rs;
	s2 = rmax/rs;
	int k;

	fseek(ff,StripOffsets[s1]+r1*nbytes*ImageWidth,SEEK_SET);
	if(s1<s2){
		fread(vec,sizeof(Ipp8u),nbytes*ImageWidth*(rs-r1),ff);
		k = rs-r1;
		for(int i=s1+1;i<s2;i++){
			fseek(ff,StripOffsets[i],SEEK_SET);
			fread(vec+k*nbytes*ImageWidth,sizeof(Ipp8u),nbytes*ImageWidth*rs,ff);
			k+=rs;
		}
		fseek(ff,StripOffsets[s2],SEEK_SET);
		fread(vec+k*nbytes*ImageWidth,sizeof(Ipp8u),nbytes*ImageWidth*(r2+1),ff);
	}else{
		fread(vec,sizeof(Ipp8u),nbytes*ImageWidth*(r2-r1+1),ff);
	}
	if(!isII){
		if(BitsPerSample == 16){
			ippsSwapBytes_16u_I((Ipp16u*)vec, (rmax-rmin+1)*ImageWidth);
		}else if(BitsPerSample == 32){
			ippsSwapBytes_32u_I((Ipp32u*)vec, (rmax-rmin+1)*ImageWidth);
		}
	}


		
	return true;
}


bool gtiff::read_int_16u(char *fileName, Ipp16u *vec){
	int nt, buff_size, iup, rem_buff_size;

	if(ff != NULL)fclose(ff);
	ff = fopen(fileName,"rb");
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	long new_file_size;
	fseek(ff, 0, SEEK_END);
	new_file_size = ftell(ff);
	
	if(new_file_size <=10){	
		sprintf(err_message,"file size (%li)",new_file_size);
		return false;
	}

	if(vect != NULL || new_file_size != file_size){
		ippsFree(vect); vect = NULL;
		file_size = new_file_size;
	}

	vect = ippsMalloc_8u(file_size);
	fseek(ff,0,SEEK_SET);
	fread(vect,sizeof(Ipp8u),file_size,ff);
	if(!fillTags()){
		ippsFree(vect); vect = NULL;
		return false;
	}
	
	nt=ImageWidth*ImageLength;

	buff_size=8*StripByteCounts[0]/BitsPerSample;
	iup = (int)StripsPerImage;
	rem_buff_size = nt%buff_size;
	

	if(rem_buff_size == 0){
		rem_buff_size = buff_size;
		if(nt>buff_size)iup--;
	}

	switch(BitsPerSample){
		case 8:
			for(int i=0;i<iup-1;i++){
				ippsConvert_8s16s((Ipp8s*)(vect+StripOffsets[i]), (Ipp16s*)(vec+i*buff_size), buff_size);
			}
			ippsConvert_8s16s((Ipp8s*)(vect+StripOffsets[iup-1]), (Ipp16s*)(vec+(iup-1)*buff_size), rem_buff_size);
			break;
		case 16:
			for(int i=0;i<iup-1;i++){
				if(!isII)ippsSwapBytes_16u_I((Ipp16u*)(vect+StripOffsets[i]), buff_size);
				ippsCopy_8u(vect+StripOffsets[i], (Ipp8u*)(vec+i*buff_size), 2*buff_size);
			}
			if(!isII)ippsSwapBytes_16u_I((Ipp16u*)(vect+StripOffsets[iup-1]), buff_size);
			ippsCopy_8u(vect+StripOffsets[iup-1], (Ipp8u*)(vec+(iup-1)*buff_size), 2*rem_buff_size);
			break;

		default:
			ippsZero_8u((Ipp8u*)vec,2*nt);
			break;
	}
	return true;
	
}

bool gtiff::read_int_8u(char *fileName, Ipp8u *vec){
	int nt, buff_size, iup, rem_buff_size;

	if(ff != NULL)fclose(ff);
	ff = fopen(fileName,"rb");
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	long new_file_size;
	fseek(ff, 0, SEEK_END);
	new_file_size = ftell(ff);
	
	if(new_file_size <=10){	
		sprintf(err_message,"file size (%li)",new_file_size);
		return false;
	}

	if(vect != NULL || new_file_size != file_size){
		ippsFree(vect); vect = NULL;
		file_size = new_file_size;
	}

	vect = ippsMalloc_8u(file_size);
	fseek(ff,0,SEEK_SET);
	fread(vect,sizeof(Ipp8u),file_size,ff);
	if(!fillTags()){
		ippsFree(vect); vect = NULL;
		return false;
	}
	
	nt=ImageWidth*ImageLength;

	buff_size=8*StripByteCounts[0]/BitsPerSample;
	iup = (int)StripsPerImage;
	rem_buff_size = nt%buff_size;
	

	if(rem_buff_size == 0){
		rem_buff_size = buff_size;
		if(nt>buff_size)iup--;
	}

	switch(BitsPerSample){
		case 8:
			for(int i=0;i<iup-1;i++){
				ippsCopy_8u(vect+StripOffsets[i], vec+i*buff_size, buff_size);
			}
			ippsCopy_8u(vect+StripOffsets[iup-1], vec+(iup-1)*buff_size, rem_buff_size);
			break;

		default:
			ippsZero_8u(vec,nt);
			break;
	}
	return true;
	
}


bool gtiff::print_short(Ipp16u a){
	*(Ipp16u*)(vect+cur_pos) = a;
	if(isII == false)ippsSwapBytes_16u_I((Ipp16u*)(vect+cur_pos),1);
	cur_pos+=2;
	return true;	
}

bool gtiff::print_long(Ipp32u a){
	*(Ipp32u*)(vect+cur_pos) = a;
	if(isII == false)ippsSwapBytes_32u_I((Ipp32u*)(vect+cur_pos),1);
	cur_pos+=4;
	return true;	
}


bool gtiff::print_ifd(Ipp16u tag, Ipp16u type, Ipp32u count, Ipp32u value){
	print_short(tag);
	print_short(type);
	print_long(count);
	print_long(value);
	return true;
}


bool gtiff::fillTemplate(){
	if(isII){
		vect[0] = 73;
		vect[1] = 73;
	}else{
		vect[0] = 77;
		vect[1] = 77;
	} 
	cur_pos = 2;
	print_short(42);
	print_long(8);

	print_short(Ipp16u(num_entries));

	StripsPerImage=1; 
	RowsPerStrip = ImageLength;
	StripByteCountsPerSample=BitsPerSample/8;

	StripOffsets=ippsMalloc_32u(StripsPerImage);
	StripByteCounts=ippsMalloc_32u(StripsPerImage);

	Ipp32u StripOffsets_min=14+12*num_entries+16;
	if(SamplesPerPixel>2){
		StripOffsets_min+=(2*(SamplesPerPixel-1));
	}
	
	StripOffsets[0]=StripOffsets_min;
	StripByteCounts[0]=StripByteCountsPerSample*ImageWidth*ImageLength;
	
	OffsetStripOffsets=StripOffsets[0];
	OffsetStripByteCounts=StripByteCounts[0];
	OffsetXResolution=14+12*num_entries;
	OffsetYResolution=OffsetXResolution+8;
	OffsetBitsPerSample = OffsetXResolution+16;
	ResolutionUnit=1;

	print_ifd(256,3,1,ImageWidth);
	print_ifd(257,3,1,ImageLength);
	if(SamplesPerPixel == 1){
		print_ifd(258,3,1,BitsPerSample);
	}else{
		print_ifd(258,3,SamplesPerPixel,OffsetBitsPerSample);
	}
	print_ifd(259,3,1,Compression);
	print_ifd(262,3,1,PhotometricInterpretation);
	print_ifd(273,4,StripsPerImage,StripOffsets[0]);
	print_ifd(277,3,1,SamplesPerPixel);
	print_ifd(278,3,1,RowsPerStrip);
	print_ifd(279,4,StripsPerImage,StripByteCounts[0]);
	print_ifd(282,5,1,OffsetXResolution);
	print_ifd(283,5,1,OffsetYResolution);
	print_ifd(296,3,1,ResolutionUnit);
	print_ifd(339,3,1,SampleFormat);
		
	print_long(0);
	
	print_long(XResolution[0]);
	print_long(XResolution[1]);
	print_long(YResolution[0]);
	print_long(YResolution[1]);
		
	if(SamplesPerPixel>2){
		for(Ipp32u i=0; i<SamplesPerPixel; i++){
			print_short((Ipp16u)BitsPerSample);
		}
	}
	return true;
}


bool gtiff::print_float(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny, unsigned int bits){
	ff = fopen(fileName,"wb");
	int nt = nx * ny;
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	int memsize;
	int nb;

	switch(bits){
		case 8:
			BitsPerSample = 8;
			break;
		case 16:
			BitsPerSample = 16;
			break;
		case 32:
			BitsPerSample = 32;
			SampleFormat=3;
			break;
		default:
			setError("wrong number of BitsPerSample");
			return false;
	}
	nb = BitsPerSample/8;
	//num_entries = 12;
	memsize = 30+12*num_entries+nb*nt;

	if(!(ImageWidth == nx && ImageLength == ny)){
		if(vect != NULL)	ippsFree(vect);	
		ImageWidth = nx;
		ImageLength = ny;
		vect = ippsMalloc_8u(memsize);
		if(vect == NULL){
			setError("can not allocate memory");
			return false;
		}
	}

	cur_pos = 0;
	fillTemplate();
	switch(bits){
		case 8:
			ippsConvert_32f8u_Sfs(vec, vect+cur_pos, nt, ippRndNear, 0);
			break;
		case 16:
			ippsConvert_32f16u_Sfs(vec, (Ipp16u*)(vect+cur_pos), nt, ippRndNear, 0);
			break;
		case 32:
			ippsCopy_32f(vec, (Ipp32f*)(vect+cur_pos),nt);
			break;
	}

	fwrite(vect,sizeof(Ipp8u),memsize,ff);
	fclose(ff);
	ff = NULL;
	return true;
}


bool gtiff::print_int_16u(char *fileName, Ipp16u *vec, unsigned int nx, unsigned int ny){
	ff = fopen(fileName,"wb");
	int nt = nx * ny;
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	int memsize;
	int nbytes;
	
	BitsPerSample = 16;
	
	nbytes = BitsPerSample/8;
	//num_entries = 12;
	memsize = 30+12*num_entries+nbytes*nt;

	if(!(ImageWidth == nx && ImageLength == ny)){
		if(vect != NULL) ippsFree(vect);	
		ImageWidth = nx;
		ImageLength = ny;
		vect = ippsMalloc_8u(memsize);
		if(vect == NULL){
			setError("can not allocate memory");
			return false;
		}
	}

	cur_pos = 0;
	fillTemplate();
		
	ippsCopy_8u((Ipp8u*)vec, vect+cur_pos, nbytes*nt);
	fwrite(vect,sizeof(Ipp8u),memsize,ff);
	fclose(ff);
	ff = NULL;
	return true;
}


bool gtiff::print_int_32u(char *fileName, Ipp32u *vec, unsigned int nx, unsigned int ny){
	ff = fopen(fileName,"wb");
	int nt = nx * ny;
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	int memsize;
	int nbytes;
	
	BitsPerSample = 32;
	
	nbytes = BitsPerSample/8;
	//num_entries = 12;
	memsize = 30+12*num_entries+nbytes*nt;

	if(!(ImageWidth == nx && ImageLength == ny)){
		if(vect != NULL) ippsFree(vect);	
		ImageWidth = nx;
		ImageLength = ny;
		vect = ippsMalloc_8u(memsize);
		if(vect == NULL){
			setError("can not allocate memory");
			return false;
		}
	}

	cur_pos = 0;
	fillTemplate();
		
	ippsCopy_8u((Ipp8u*)vec, vect+cur_pos, nbytes*nt);
	fwrite(vect,sizeof(Ipp8u),memsize,ff);
	fclose(ff);
	ff = NULL;
	return true;
}


bool gtiff::print_int_8u(char *fileName, Ipp8u *vec, unsigned int nx, unsigned int ny){
	ff = fopen(fileName,"wb");
	int nt = nx * ny;
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	int memsize;
	int nbytes;
	
	BitsPerSample = 8;
	
	nbytes = BitsPerSample/8;
	//num_entries = 12;
	memsize = 30+12*num_entries+nbytes*nt;

	if(!(ImageWidth == nx && ImageLength == ny)){
		if(vect != NULL) ippsFree(vect);	
		ImageWidth = nx;
		ImageLength = ny;
		vect = ippsMalloc_8u(memsize);
		if(vect == NULL){
			setError("can not allocate memory");
			return false;
		}
	}

	cur_pos = 0;
	fillTemplate();
		
	ippsCopy_8u(vec, vect+cur_pos, nbytes*nt);
	fwrite(vect,sizeof(Ipp8u),memsize,ff);
	fclose(ff);
	ff = NULL;
	return true;
}


bool gtiff::print_rgb_8u(char *fileName, Ipp8u *vec, unsigned int nx, unsigned int ny){
	ff = fopen(fileName,"wb");
	int nt = nx * ny;
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}

	int memsize;
	int nbytes;
	
	PhotometricInterpretation=2;
	BitsPerSample = 8;
	
	SamplesPerPixel = 3;
	nbytes = SamplesPerPixel*BitsPerSample/8;
	//num_entries = 13;
	memsize = 36+12*num_entries+nbytes*nt;

	if(!(ImageWidth == nx && ImageLength == ny)){
		if(vect != NULL) ippsFree(vect);	
		ImageWidth = nx;
		ImageLength = ny;
		vect = ippsMalloc_8u(memsize);
		if(vect == NULL){
			setError("can not allocate memory");
			return false;
		}
	}

	cur_pos = 0;
	fillTemplate();
		
	ippsCopy_8u(vec, vect+cur_pos, nbytes*nt);
	fwrite(vect,sizeof(Ipp8u),memsize,ff);
	fclose(ff);
	ff = NULL;
	return true;
}



bool gtiff::print_float_16u(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny, float vmin, float vmax){
	int nt = nx * ny;
	float ratio = 65535.0f/(vmax-vmin);
	ippsSubC_32f_I(vmin,vec,nt);
	ippsMulC_32f_I(ratio,vec,nt);
	return print_float_16u(fileName, vec, nx, ny);
}


bool gtiff::print_float_8u(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny, float vmin, float vmax){
	int nt = nx * ny;
	float ratio = 255.0f/(vmax-vmin);
	ippsSubC_32f_I(vmin,vec,nt);
	ippsMulC_32f_I(ratio,vec,nt);
	return print_float_8u(fileName, vec, nx, ny);
}

bool gtiff::print_float_32f(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny, float vmin, float vmax){
	int nt = nx * ny;
	ippsThreshold_LTValGTVal_32f_I(vec,nt,vmin,vmin,vmax,vmax);
	return print_float_32f(fileName, vec, nx, ny);
}

bool gtiff::print_float_8u(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny){
	return print_float(fileName, vec, nx, ny, 8);
}

bool gtiff::print_float_16u(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny){
	return print_float(fileName, vec, nx, ny, 16);
}


bool gtiff::print_float_32f(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny){
	return print_float(fileName, vec, nx, ny, 32);
}


