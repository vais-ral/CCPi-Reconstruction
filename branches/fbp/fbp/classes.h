class gpu_info{
public:
	gpu_info(){};
	unsigned int sharedMemPerBlock, sharedMemBanks;
	unsigned int regsPerBlock;
	unsigned int warpSize;
	unsigned int maxResidentWarps, maxResidentThreads, maxResidentBlocks;
	unsigned int maxThreadsPerBlock;
	unsigned int major, minor;
	unsigned int multiProcessorCount;
	unsigned int nxi, nyi, nxo, nyo;
	unsigned int mxi, myi, mxo, myo;
	unsigned int wo, ho;
	float xmin, xmax, ymin, ymax, angle;
	float origWidth, origHeight;
	float newWidth, newHeight, newX, newY;
	float imageCentre, outputPixelSize, outputBlockDiameter, outputROIRadius;
	unsigned int x_min, x_max;
	float roiR, roiA;
	unsigned int vertH;
	unsigned int fpower, flength, origw;
	float rotAngleStep;
	unsigned int ux, uy, nuc;
	float blockRadius;
	unsigned int pol_r, pol_a;
	

};

class time_stamp{
public:
	time_stamp(){};
	char cname[100];
	void read_time(Ipp8u *q);
	float sec;
	unsigned int minutes, hours, days, months, years;
};


class dir_entry{
public:
	void print_dir_entry(FILE *fp);
	char name[32];
	char comm[500];
	//char time_mod[100], time_creation[100];
	wchar_t lname[32];
	int size_area;
	Ipp32s dirID_left, dirID_right, dirID_root;
	Ipp8u type, node_colour;
	long long uID[2];
	int flag;
	Ipp32s secID, stream_size;
	int vtype;
	int esize;
	time_stamp t_creation, t_modification;
	bool pnote, psize, ptime_mod, ptime_creation;
	dir_entry();
};

class xXrm2tif{
public:
	xXrm2tif();
	bool readXml(Ipp8u *vf, int s1, int s2);
	char folder[MAX_FOLDER], logFile[MAX_FOLDER], tagFile[MAX_FOLDER];
	char prefix[MAX_FILE], suffix[MAX_FILE], extension[MAX_EXT];
	char flatName[MAX_FOLDER], darkName[MAX_FOLDER], inputFile[MAX_FOLDER];
	int NOD, fileFirst, fileStep;
	int imageFirst, imageLast, imageStep;
	int imagesPerFile;
	int memorySizeMin, memorySizeMax;
	float memorySizeFMax;
	int output_type;
	int rowFirst, rowLast;
	int width, height, NOI, bitsPerData, bitsPerRef;
	float pixelSize, d2a, s2a;
	bool isRef;
};

class xAllFDK{
public:
	xAllFDK(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	char inputFolder[MAX_FOLDER], logFile[MAX_FOLDER];
	char inputPrefix[MAX_FILE], inputSuffix[MAX_FILE], inputExtension[MAX_EXT];
	char outputFolder[MAX_FOLDER];
	char outputPrefix[MAX_FILE], outputSuffix[MAX_FILE], outputExtension[MAX_EXT];
	char flatName[MAX_FOLDER], darkName[MAX_FOLDER];
	char xShiftsFile[MAX_FOLDER], yShiftsFile[MAX_FOLDER];
	char filterName[30];

	int numberOfCircles;

	int inputNOD, outputNOD;
	int outputBits;
	int inputFileFirst, inputFileStep;
	int inputImageFirst, inputImageLast, inputImageStep;
	unsigned int inputImagesPerFile;
	int flatType, darkType, inputRestrictions;
	
	float flatValue, darkValue, flatMin, darkMin, flatMax, darkMax;
	float inputMin, inputMax, outputMin, outputMax;
	int rowFirst;
	int ringArtefactsType;
	float ringArtefactsParamR, ringArtefactsParamN;
	float sliceFirst, sliceLast, sliceStep;
	int outputIndexFirst, outputIndexStep;
	float pixelSize, angleStep;
	int xShiftsType, yShiftsType, xShiftsAbs, yShiftsAbs;
	
	float sourceXShift, sourceYShift, originalImageHeight, axisXShift;
	float sourceToObject, detectorToObject;
	float rotationAngle;
	
	float filterParam;
	float ROIXmin, ROIXmax, ROIYmin, ROIYmax;
	int outputWidthType, outputWidth;
	int clockwiseRotation;
};


class xBeamlineUser{
public:
	xBeamlineUser(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	int type, year, month, day;
	char beamlineName[MAX_FILE], visitNumber[MAX_FILE];
};

class xRaw{
public:
	xRaw(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	int type, bits, offset, byteOrder, gap, xlen, ylen, zlen;
};

  
class xInputData{
public:
	xInputData();
	~xInputData();
	xRaw *raw;
	bool readXml(Ipp8u *vf, int s1, int s2);
	char folder[MAX_FOLDER], prefix[MAX_FILE], suffix[MAX_FILE], extension[MAX_EXT];
	int NOD, fileFirst, fileLast, fileStep, imageFirst, imageLast, imageStep;
	int firstImageIndex, imagesPerFile;
	int memorySizeMin, memorySizeMax;
	int type, restrictions, orientation, shape;
	float memorySizeFMax, valueMin, valueMax;
	float pixelParam;
};


class xFDField{
public:
	
	bool readXml(Ipp8u *vf, int s1, int s2);
	char fileBefore[MAX_FOLDER], fileAfter[MAX_FOLDER], fileProfile[MAX_FOLDER];
	int type, typeProfile;
	float valueBefore, valueAfter;
	xFDField(){
		fileBefore[0] = '\0';
		fileAfter[0] = '\0';
		fileProfile[0] = '\0';
	};
};

class xFlatDarkFields{
public:
	xFlatDarkFields();
	~xFlatDarkFields();
	xFDField *flatField, *darkField;
	bool readXml(Ipp8u *vf, int s1, int s2);
};

class xHighPeaks{
public:
	xHighPeaks(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	int type, numberPixels;
	float jump;
};

class xRingArtefacts{
public:
	xRingArtefacts(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	int type;
	float parameterN, parameterR;
	int num_series;
};

class xIntensity{
public:
	xIntensity(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	int type, columnLeft, columnRight, zeroLeft, zeroRight;
};


class xPreprocessing{
public:
	xPreprocessing();
	~xPreprocessing();
	bool readXml(Ipp8u *vf, int s1, int s2);
	xHighPeaks *highPeaksBefore, *highPeaksCols, *highPeaksRows;
	xRingArtefacts *ringArtefacts;
	xIntensity *intensity;
};



class xTransform{
public:
	xTransform(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	char missedProjections[MAX_FOLDER];
	int missedProjectionsType, rotationAngleType, rotationAngleEndPoints, scaleType, extrapolationType;
	float rotationAngle, reCentreAngle, reCentreRadius;
	unsigned int cropTop, cropBottom, cropLeft, cropRight;
	unsigned int scaleWidth, scaleHeight, extrapolationWidth, extrapolationPixels;
	int interpolation;

};

class xFilter{
public:
	xFilter(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	int type;
	int name, windowName;
	float bandwidth;
	float pixelSize;
	int norm;
};

class xTilt{
public:
	xTilt(){};
	int type;
	float z_tilt, x_tilt;
};

class xCoordinateSystem{
public:
	xCoordinateSystem(){};
};

class xRadius{
public:
	xRadius(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	int type;
	float percent, pixel;
};

class xCircles{
public:
	xCircles();
	~xCircles();
	bool readXml(Ipp8u *vf, int s1, int s2);
	xRadius *radiusMin, *radiusMax, *radiusStep;
};

class xRoi{
public:
	xRoi(){};
	bool readXml(Ipp8u *vf, int s1, int s2);
	float xmin, xmax, ymin, ymax;
	int type, outputWidth, outputWidthType;
	float angle; 
	
};

class xBackprojection{
public:
	xBackprojection();
	~xBackprojection();
	bool readXml(Ipp8u *vf, int s1, int s2);
	xFilter *filter;
	float imageCentre;
	xTilt *tilt;
	xCoordinateSystem *coordinateSystem;
	xCircles *circles;
	xRoi *roi;
	int pc_interpolation;
};


class xOutputData{
public:
	xOutputData();
	bool readXml(Ipp8u *vf, int s1, int s2);
	char folder[MAX_FOLDER], prefix[MAX_FILE], suffix[MAX_FILE], extension[MAX_EXT];
	int type, NOD, fileFirst, fileStep, bits, bitsType, restrictions;
	float valueMin, valueMax;
	int shape, state;
};


class xFBP{
public:
	xFBP();
	~xFBP();
	bool readXml(Ipp8u *vf, int s1, int s2);
	char defaultXml[MAX_FOLDER], logFile[MAX_FOLDER];
	int GPUDeviceNumber;
	xBeamlineUser *beamlineUser;
	xInputData *inputData;
	xFlatDarkFields *flatDarkFields;
	xPreprocessing *preprocessing;
	xTransform *transform;
	xBackprojection *backprojection;
	xOutputData *outputData;
};

class xHMset{
public:
	xHMset();
	~xHMset();
	void **pointers;
	int child_name[MAX_CHILD];
	int lenf;
	int nchild;
	
	char errorString[MAX_ERROR];
	bool readXml(char *inputFile);
	bool outState;
	xFBP *fbp;
	xXrm2tif *x2t;
	xAllFDK *af;
	Ipp8u *vf;
};

class xFBP_gpu{
public:
	xFBP_gpu(){};
	unsigned int nx, ny, nxo, nyo;
	unsigned int chunkLeft, chunkWidth, numChunks;
	unsigned int blockWidth, blockHeight, vH;
	float xc;
	//size_t;
};

class xData{
public:
	xData();
	~xData();
	void freeMemory();
	Ipp32f *matd1, *matd2, *matf1, *matf2;
	Ipp32f *profd, *proff, *veci;
	Ipp8u *vecio;
	Ipp32f *vecd, *vecf;
	Ipp32f *vect1, *vect2;
	Ipp32f *veco, *vecp1, *vecp2;
	Ipp32f *vec32;
	Ipp64f *vec64, *vec_res;
	Ipp64f *mata;
	Ipp16u *vmiss;
	Ipp32f *vtb, *vta;
	Ipp32f *mapx, *mapy;
	Ipp64f *pp_a, *pp_b, *pp_c, *pp_f, *pp_t;
	Ipp32f *vto, *vfilter;
	Ipp32f *vecArcTan;
	Ipp32f *bpx, *bpy;
	Ipp32f *vecCos;
	Ipp64f *vecRS;
	Ipp32f *shift_x, *shift_y;
	Ipp32f *vecXY, *vecbXY, *vecXY_block;

	Ipp32f *vecSinCos;
	Ipp32f *vta_warp;

	Ipp32f *vecRA, *vecAT;
	Ipp32f *vecnCos;
	Ipp32fc *veccF, *veccI;
	Ipp32f *vecX, *vecY;
	Ipp32f *vecSin, *vecCS;
	Ipp32f *vecti, *vecto, *vectr;
	Ipp32f *vecPol;
};


class xTime{
public:
	xTime(){};
	clock_t  global_start, global_end, local_end;
	clock_t local_start[MAX_TAG_COUNT];
};


class xtiff{
public:
	xtiff();
	~xtiff();
	char err_message[50], ctemp[50];
	Ipp8u *vect;
	long file_size;
	FILE *ff;
	Ipp32u cur_pos;
	
	bool isII;
	Ipp32u BitsPerSample, Compression, PhotometricInterpretation;
	Ipp32u ImageWidth, ImageLength, RowsPerStrip;
	Ipp32u ResolutionUnit;
	Ipp32u XResolution[2], YResolution[2];
	Ipp32u OffsetXResolution, OffsetYResolution, OffsetStripByteCounts;
	Ipp32u OffsetStripOffsets;
	Ipp32u *StripOffsets;
	Ipp32u *StripByteCounts;
	Ipp32u SamplesPerPixel, OffsetBitsPerSample;
	Ipp32u StripsPerImage;
	Ipp32u StripByteCountsPerSample;
	Ipp32u SampleFormat;

	Ipp32u num_entries;
	
	bool getInfo(Ipp8u *vecio, long file_size, unsigned int *nx, unsigned int *ny, unsigned int *nb);
	bool read_float(Ipp8u *vecio, long ifile_size, Ipp32f *vec);
	long get_file_size(unsigned int nx, unsigned int ny, unsigned int nb);

	void setError(char *err_string);
	bool fillTags();

	bool read_short(Ipp16u *);
	bool read_long(Ipp32u *);
	bool fread_short(Ipp16u *);
	bool fread_long(Ipp32u *);
	bool find_next_ifd(Ipp32u *off);
	bool ffind_next_ifd(Ipp32u *off);
	bool read_entry();
	bool fread_entry();

	bool print_int_32u(char *fileName, Ipp32u *vec, unsigned int nx, unsigned int ny);
	bool print_int_16u(char *fileName, Ipp16u *vec, unsigned int nx, unsigned int ny);
	bool print_int_8u(char *fileName, Ipp8u *vec, unsigned int nx, unsigned int ny);
	bool print_rgb_8u(char *fileName, Ipp8u *vec, unsigned int nx, unsigned int ny);

	bool print_float(Ipp8u *vecio, Ipp32f *vec, unsigned int nx, unsigned int ny, unsigned int bits);

	bool print_ifd(Ipp16u tag, Ipp16u type, Ipp32u count, Ipp32u value);
	bool print_short(Ipp16u a);
	bool print_long(Ipp32u a);
	bool fillTemplate();

	
	
};


class gtiff{
public:
	gtiff();
	~gtiff();
	char err_message[50], ctemp[50];
	Ipp8u *vect;
	long file_size;
	FILE *ff;
	Ipp32u cur_pos;
	
	void setError(char *err_string);
	bool fillTags();
	bool ffillTags();
	bool getInfo(char *fileName, unsigned int *nx, unsigned int *ny, unsigned int *nb);

	bool read_short(Ipp16u *);
	bool read_long(Ipp32u *);
	bool fread_short(Ipp16u *);
	bool fread_long(Ipp32u *);
	bool find_next_ifd(Ipp32u *off);
	bool ffind_next_ifd(Ipp32u *off);
	bool read_entry();
	bool fread_entry();
	bool read_float(char *fileName, Ipp32f *vec);
	bool read_int_16u(char *fileName, Ipp16u *vec);
	bool read_int_8u(char *fileName, Ipp8u *vec);
	bool read_chunk(char *fileName, Ipp8u *vec, int rmin, int rmax);

	bool print_int_32u(char *fileName, Ipp32u *vec, unsigned int nx, unsigned int ny);
	bool print_int_16u(char *fileName, Ipp16u *vec, unsigned int nx, unsigned int ny);
	bool print_int_8u(char *fileName, Ipp8u *vec, unsigned int nx, unsigned int ny);
	bool print_rgb_8u(char *fileName, Ipp8u *vec, unsigned int nx, unsigned int ny);

	bool print_float(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny, unsigned int bits);
	bool print_float_8u(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny, float vmin, float vmax);
	bool print_float_8u(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny);
	bool print_float_16u(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny, float vmin, float vmax);
	bool print_float_16u(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny);
	bool print_float_32f(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny, float vmin, float vmax);
	bool print_float_32f(char *fileName, Ipp32f *vec, unsigned int nx, unsigned int ny);

	bool print_ifd(Ipp16u tag, Ipp16u type, Ipp32u count, Ipp32u value);
	bool print_short(Ipp16u a);
	bool print_long(Ipp32u a);
	bool fillTemplate();

	bool isII;
	Ipp32u BitsPerSample, Compression, PhotometricInterpretation;
	Ipp32u ImageWidth, ImageLength, RowsPerStrip;
	Ipp32u ResolutionUnit;
	Ipp32u XResolution[2], YResolution[2];
	Ipp32u OffsetXResolution, OffsetYResolution, OffsetStripByteCounts;
	Ipp32u OffsetStripOffsets;
	Ipp32u *StripOffsets;
	Ipp32u *StripByteCounts;
	Ipp32u SamplesPerPixel, OffsetBitsPerSample;
	Ipp32u StripsPerImage;
	Ipp32u StripByteCountsPerSample;
	Ipp32u SampleFormat;

	Ipp32u num_entries;
	
};


class MatPar{
public:
	MatPar();
	~MatPar();
	bool formCudaMatrixXY();
	float rmin, rmax;
	int na, nr;
	int ImageWidth, ImageHeight;
	//float xc, yc;
	float y_min, y_max;
	int y_mini, y_maxi;
	float SourceToDetector, PixelSize;
	float rowf;
	float ConeCentreShiftX, ConeCentreShiftY;
	float AxisShift;
	int RestrMin, RestrMax, RowFirst;
	Ipp32f *m_x, *m_y;
	Ipp32f *ivol, *rec;
	Ipp32fc *vsect, *fsect;
	//float prmax, prmin;
};

class OtherParam{
public:
	int l_a, l_r, l_s;
	int w_a, w_r, w_s;
	int num_images, images_to_process;
	int eWidth, eWidthP;
	Ipp32f *vecFilterRowR;
	OtherParam(){
		vecFilterRowR = NULL;
	};

	~OtherParam(){
		if(vecFilterRowR != NULL){ippsFree(vecFilterRowR); vecFilterRowR = NULL;}
	};
};


class allData{
public:
	allData(xHMset *hmset);
	~allData();
	
	char ifile[MAX_FOLDER], ofile[MAX_FOLDER];
	char itemplate[MAX_FOLDER], otemplate[MAX_FOLDER];
	char tagName[MAX_TAG_LENGTH * MAX_TAG_COUNT], message[MAX_MESSAGE], cspace[MAX_FILE];

	xHMset *hms;
	xData *data;
	xTime *time;
	gpu_info *gi;
	xFBP_gpu *fg;
	FILE* fp;

	unsigned int row_first, row_last, row_step, row;
	int nz, nt, nr;
	unsigned int nx, ny, nb;
	int nx_out, ny_out, nb_out;
	unsigned int dark_nx, dark_ny, flat_nx, flat_ny;
	unsigned int profd_nx, profd_ny, proff_nx, proff_ny;
	int miss_proj;
	int tb_nx, tb_ny;
	int ta_nx, ta_ny, ta_ny_13;
	int to_nx, to_ny;
	int ixc, dxc, bpw, sxc;
	int tagNumber;
	int ny_warp;
	unsigned int temp_nx, temp_ny, temp_nb;
	float new_xc;
	
	float output_x_min, output_x_max, output_y_min, output_y_max;
	float roiR, roiA, roiX, roiY, roiRR, roiAA, scale_factor; 
	bool isFirstSlice, isReCentre;

	float roiRadius;
};


class PolarToCart{
public:
	PolarToCart();
	~PolarToCart();
	bool formCudaMatrix();
	bool formCudaMatrix(allData *aD);
	int ImageWidth;
	int OutputWidth, OutputHeight;
	int NumberOfCircles;
	int num_images;
	float xmin, xmax, ymin, ymax;
	float rmin, rmax;
	float RotAngle;
	int ailen, aimin;
	int l_a, l_r;
	int w_a, w_r;
	int nan, nrn;
	Ipp32f *vec_x, *vec_y, *vec_a, *vec_r;
	Ipp32s ntb;
	Ipp32s *i_a, *i_r;
	Ipp32f *pc_r, *pc_a;
	Ipp32s *pc_pos, *pc_size;
	Ipp32f *pc_res;
	Ipp32s pc_len;
	int pc_nc;
};



class ole2cd{
public:
	ole2cd();
	~ole2cd();
	FILE *ff, *fo;
	char cdfID[9], UID[17];
	Ipp8u q[128];
	Ipp16u revision_number, version_number;
	Ipp32u num_sect_dir, num_sec_for_SAT, SecID_first_dir, min_size_std_stream;
	Ipp32u size_sector, size_short_sector, num_SAT, SecID_first_SSAT, num_SSAT;
	Ipp32u SecID_first_MSAT, num_MSAT, tot_dir;
	
	int root_ind;
	Ipp32s start_ss;

	dir_entry **p_dir;
	int *num_child, *child_ind, *parent;
	int *vvv, *emb, *emb_ind;
	int *sectors;
	Ipp32u *MSAT, *SSAT;
	Ipp32s *SAT;

	bool byte_order;
	int pow2(int s);
	bool read_short(int *s);
	bool read_bool(int *s);
	bool read_long(int *s);
	bool read_float(float *f);
	
	bool read_types(allData *aD);
	void print_xml_file();
	bool print_xml(allData *aD, char *sw);
	bool read_header();
	bool sort_dir();
	void read_dir_entry(dir_entry *p);

	bool fill_MainParam(allData *aD);
	bool print_MainParam(allData *aD, char *finfo);

	bool print_data(allData *aD, char *fraw);
	bool print_ref_data(allData *aD);
};


