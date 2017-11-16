/* ----------------------------------------------------------------------
 * xtiff.cpp()
 * ----------------------------------------------------------------------
 * Read and write a small set of tiff file formats.
 *
 * TODO: Merge xtiff and gtiff classes. Interface to libtiff.
 * Valeriy: gtiff will be removed later. Now it is used for "cone" functions only.
 * ----------------------------------------------------------------------
 */
#include "defs.h"

xtiff::xtiff(){
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

// Read the TIFF Image File Directory and parse the tags in there
bool xtiff::find_next_ifd(Ipp32u *off){
	Ipp16u num_of_dir;
	read_short(&num_of_dir);
	num_entries += Ipp32u(num_of_dir);

	for(Ipp16u i=0; i<num_of_dir; i++)
		read_entry();

	read_long(off);
	return true;
}


bool xtiff::ffind_next_ifd(Ipp32u *off){
	Ipp16u num_of_dir;
	fread_short(&num_of_dir);
	num_entries += Ipp32u(num_of_dir);

	for(Ipp16u i=0;i<num_of_dir;i++)
		fread_entry();

	fread_long(off);
	return true;
}


void xtiff::setError(char *err_string){
	sprintf(err_message,err_string);
	printf("Error: %s.\n",err_message);
}

bool xtiff::getInfo(Ipp8u *vecio, long ifile_size, unsigned int *nx, unsigned int *ny, unsigned int *nb){
	file_size = ifile_size;
	vect = vecio;
	if(!fillTags()){
		return false;
	}
	
	*nx = ImageWidth;
	*ny = ImageLength;
	*nb = BitsPerSample;

	return true;
}

bool xtiff::fread_short(Ipp16u *a){
	fread(a,sizeof(Ipp16u),1,ff);
	//*a = *((Ipp16u*)(vect+cur_pos));
	if(isII == false)ippsSwapBytes_16u_I(a,1);
	cur_pos+=2;
	return true;
}

bool xtiff::read_short(Ipp16u *a){
	*a = *((Ipp16u*)(vect+cur_pos));
	if(isII == false)ippsSwapBytes_16u_I(a,1);
	cur_pos+=2;
	return true;
}

bool xtiff::fread_long(Ipp32u *a){
	fread(a,sizeof(Ipp32u),1,ff);
//	*a = *((Ipp32u*)(vect+cur_pos));
	if(isII == false)ippsSwapBytes_32u_I(a,1);
	cur_pos+=4;
	return true;
}


bool xtiff::read_long(Ipp32u *a){
	*a = *((Ipp32u*)(vect+cur_pos));
	if(isII == false)ippsSwapBytes_32u_I(a,1);
	cur_pos+=4;
	return true;
}


bool xtiff::fillTags(){
	
  // First two bytes of a TIFF file represent the characters II or MM:
  // II for little-endian, MM for big-endian.

	if(vect[0] != vect[1]){
		sprintf(ctemp,"this is not TIFF, first two bytes are %i, %i",vect[0],vect[1]);
		setError(ctemp);
		return false;
	}
	if(vect[0] == 73){
		isII = true;	// II
	}else  if(vect[1] == 77){
		isII = false;	// MM
	}else{
		sprintf(ctemp,"this is not TIFF, first two byte are %i, %i",vect[0],vect[1]);
		setError(ctemp);
		return false;
	}
	cur_pos = 2;

	unsigned short sh;
	read_short(&sh);

	// Next two bytes represents the number 42.
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

bool xtiff::read_entry(){
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
			  BitsPerSample=value/65536; // GWL: Why 65536? 
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


bool xtiff::fread_entry(){
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


xtiff::~xtiff(){
	if(ff != NULL) fclose(ff);
	if(StripOffsets != NULL){
		ippsFree(StripOffsets); StripOffsets = NULL;
	}
	if(StripByteCounts != NULL){
		ippsFree(StripByteCounts); StripByteCounts = NULL;
	}
	/*if(vect != NULL){
		ippsFree(vect); vect = NULL;
	}*/
}


bool xtiff::read_float(Ipp8u *vecio, long ifile_size, Ipp32f *vec){
	int nt, buff_size, iup, rem_buff_size;
	vect = vecio;
	file_size = ifile_size;
	
	if(!fillTags()){
		return false;
	}
	
	nt=ImageWidth*ImageLength;

	buff_size=8*StripByteCounts[0]/BitsPerSample;
	iup = (int)StripsPerImage;
	if(iup == 1) buff_size = 8*this->OffsetStripByteCounts/BitsPerSample;
	rem_buff_size = nt%buff_size;

	if(rem_buff_size == 0){
		iup++;
		//rem_buff_size = buff_size;
		//if(nt>buff_size)iup--;
	}

	switch(BitsPerSample){
		case 8:
			for(int i=0;i<iup-1;i++){
				ippsConvert_8u32f(vect+StripOffsets[i], vec+i*buff_size, buff_size);
			}
			if(rem_buff_size > 0){
				ippsConvert_8u32f(vect+StripOffsets[iup-1], vec+(iup-1)*buff_size, rem_buff_size);
			}
			break;
		case 16:
			for(int i=0;i<iup-1;i++){
				if(!isII)ippsSwapBytes_16u_I((Ipp16u*)(vect+StripOffsets[i]), buff_size);
				ippsConvert_16u32f((Ipp16u*)(vect+StripOffsets[i]), vec+i*buff_size, buff_size);
			}
			if(rem_buff_size > 0){
				if(!isII)ippsSwapBytes_16u_I((Ipp16u*)(vect+StripOffsets[iup-1]), buff_size);
				ippsConvert_16u32f((Ipp16u*)(vect+StripOffsets[iup-1]), vec+(iup-1)*buff_size, rem_buff_size);
			}
			break;
		case 32:
			for(int i=0;i<iup-1;i++){
				if(!isII)ippsSwapBytes_32u_I((Ipp32u*)(vect+StripOffsets[i]), buff_size);
				ippsCopy_32f((Ipp32f*)(vect+StripOffsets[i]), vec+i*buff_size, buff_size);
			}
			if(rem_buff_size > 0){
				if(!isII)ippsSwapBytes_32u_I((Ipp32u*)(vect+StripOffsets[iup-1]), buff_size);
				ippsCopy_32f((Ipp32f*)(vect+StripOffsets[iup-1]), vec+(iup-1)*buff_size, rem_buff_size);
			}
			break;

		default:
			ippsZero_32f(vec,nt);
			break;
	}
	
	return true;
}

bool xtiff::print_short(Ipp16u a){
	*(Ipp16u*)(vect+cur_pos) = a;
	if(isII == false)ippsSwapBytes_16u_I((Ipp16u*)(vect+cur_pos),1);
	cur_pos+=2;
	return true;	
}

bool xtiff::print_long(Ipp32u a){
	*(Ipp32u*)(vect+cur_pos) = a;
	if(isII == false)ippsSwapBytes_32u_I((Ipp32u*)(vect+cur_pos),1);
	cur_pos+=4;
	return true;	
}


bool xtiff::print_ifd(Ipp16u tag, Ipp16u type, Ipp32u count, Ipp32u value){
	print_short(tag);
	print_short(type);
	print_long(count);
	print_long(value);
	return true;
}


bool xtiff::fillTemplate(){
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
		for(unsigned int i=0;i<SamplesPerPixel;i++){
			print_short(Ipp16u(BitsPerSample));
		}
	}
	return true;
}


bool xtiff::print_float(Ipp8u *vecio, Ipp32f *vec, unsigned int nx, unsigned int ny, unsigned int bits){
	int nt;
	nt = nx * ny;
	if(nt <= 0){
		setError("wrong size");
		return false;
	}

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
	ImageWidth = nx;
	ImageLength = ny;

	vect = vecio;

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

	return true;
}


bool xtiff::print_int_16u(char *fileName, Ipp16u *vec, unsigned int nx, unsigned int ny){
	unsigned int nt;
	ff = fopen(fileName,"wb");
	
	if(ff == NULL){
		setError("can not open the file");
		return false;
	}
	nt = nx * ny;

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


bool xtiff::print_int_32u(char *fileName, Ipp32u *vec, unsigned int nx, unsigned int ny){
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


bool xtiff::print_int_8u(char *fileName, Ipp8u *vec, unsigned int nx, unsigned int ny){
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


bool xtiff::print_rgb_8u(char *fileName, Ipp8u *vec, unsigned int nx, unsigned int ny){
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



long xtiff::get_file_size(unsigned int nx, unsigned int ny, unsigned int nb){
	return 30+12*num_entries+(nb/8)*nx*ny;
}
