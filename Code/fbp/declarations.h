bool allocate_memory(allData *aD);

#if defined(_WIN32) || defined(_WIN64)
char *basename( char *path );
#endif

bool bp_all(allData *aD);

bool check_crop(allData *aD);
bool check_is_null(allData *aD, Ipp32f *v, const char *name, const char *comm);
bool check_rows(allData *aD);
bool check_roi(allData *aD);
bool check_roi_size(allData *aD);
bool closeAllTags(allData *aD);
bool closeAllTags(allData *aD, bool tres);
bool createFolder(allData *aD, char *fname);
bool createFolderCycle(allData *aD, char *fname);
int cudaRound(int v, int d);

bool fbp_gpu(allData *aD);
bool ff_correction_af(allData *aD);
bool ff_correction(allData *aD);


bool fillMatPar(MatPar *matp, PolarToCart *PC, allData *aD, OtherParam *other, float row);
bool fillPC(PolarToCart *PC, allData *aD, OtherParam *other);


bool find_missed_projections(allData *aD);
bool findEndTag(Ipp8u *vf, int poss, int pose, int *pos1, int *pos2, char *tagName);
bool findNextTag(Ipp8u *vf, int start_pos, int end_pos, int *tag_left_pos, int *tag_right_pos, char *tagName, bool *isEmpty);
bool find_power(int origw, int *next, int *extpower);
bool findSlices(int *SliceFirst, int *SliceLast, int *SliceStep, int RestrFirst, int RestrLast);
bool form_files_names(allData *aD);
bool form_files_names(allData *aD, int row);
bool form_filter(allData *aD);
bool find_fbp_param_20(allData *aD);

void formFilterRL2(Ipp32f *vec,Ipp32f omega, Ipp32f ds, int nx);
bool formRowFilter(allData *aD, OtherParam *other);

bool getValue(Ipp8u *vf, int ch2, int ch3, char *value);
bool get_gpu_info(allData *aD);

bool high_peaks_before(allData *aD);

bool ifFolderExists(allData *aD, char *fname);
void init_timestamp();
void int_ratio(unsigned int m, unsigned int d, unsigned int *n);
bool intensity_norm(allData *aD);

bool main_allfdk(allData *aD);
bool main_fbp(allData *aD);
bool main_xrm2tif(allData *aD);
bool motion_correction_af(allData *aD);

bool nom(int a, int b, int *c);


bool preprocess_input_data(allData *aD);

void printusage(int argc, char **argv);

void print_all(dir_entry *p, ole2cd *od, FILE *fp, Ipp8u *q);
void print_bottom(dir_entry *p, FILE *fp);
bool print_global_time_stamp(allData *aD);
bool print_local_time_stamp(allData *aD);
bool print_output_data(allData *aD);
void print_sh(FILE *fp, Ipp8u* q, int m, int n, int pe, int vt, char * name);

void print_top(dir_entry *p, FILE *fp);

bool printError(allData *aD);
bool printError(allData *aD, char *text);
bool printInfo(allData *aD);
bool printInfo(allData *aD, char *text);
bool printSettingsXml_FBP(allData *aD);
bool printTag(allData *aD, char *tag, int i, char *comm);
bool printTag(allData *aD, char *tag, int i);
bool printTag(allData *aD, char *tag, long i, char *comm);
bool printTag(allData *aD, char *tag, long i);
bool printTag(allData *aD, char *tag, unsigned int i, char *comm);
bool printTag(allData *aD, char *tag, unsigned int i);
bool printTag(allData *aD, char *tag, float i, char *comm);
bool printTag(allData *aD, char *tag, float i);
bool printTag(allData *aD, char *tag, char *text, char *comm);
bool printTag(allData *aD, char *tag, char *text);
bool printTagEnd(allData *aD);
bool printTagEnd2(allData *aD);
bool printTagStart(allData *aD, char *tag);
bool printTagStart(allData *aD, char *tag, char *comm);
bool printTagStart2(allData *aD, char *tag);
bool printTagStart2(allData *aD, char *tag, char *comm);
bool printWarning(allData *aD);
bool printWarning(allData *aD, char *text);

bool read_fdf_data(allData *aD);
bool read_field(allData *aD, int df);
bool read_input_data(allData *aD);
bool read_input_data_af(allData *aD);
bool read_one_slice(allData *aD);
bool read_one_slice_fs(allData *aD);
bool ring_artefacts_af(allData *aD);
bool ring_artefacts_removal(allData *aD);
bool ring_artefacts_removal_fs(allData *aD);

bool read_vector_xml(allData *aD, char *inputFile, Ipp32f *vec);

bool readInputVolume(MatPar *matp, PolarToCart *PC, allData *aD, OtherParam *other);
int removeContiguousSpace(Ipp8u *str, int *len, int start_pos);
int removeContiguousSpaceRev(Ipp8u *str, int *len, int start_pos);
bool replaceSlash(char *value);

bool scale_data(allData *aD);
bool set_local_time_stamp(allData *aD);
Ipp32f sinc(Ipp32f s);

void timestamp(const char *stampmsg, const int vlevel);
void time_to_char(time_stamp *ts, char *c);
bool toLower(char *value);
bool toUpper(char *value);
bool transform_sinogram(allData *aD);
bool transform_sinogram_fs(allData *aD);
bool transform_sinogram_sol_fs(allData *aD);
bool transform_sinogram_trans_fs(allData *aD);

bool cudaReconstructFDK(MatPar *matp, PolarToCart *PC, allData *aD, OtherParam *other);
bool fbp_axial_cuda(allData *aD);
bool fbp_cuda(allData *aD);

bool fbp_cuda_20(allData *aD);
bool fbp_cuda_cpu2(allData *aD);

bool find_params(allData *aD);
