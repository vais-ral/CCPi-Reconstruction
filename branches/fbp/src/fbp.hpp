
#ifndef CCPI_FBP_FILTERS
#define CCPI_FBP_FILTERS

enum filter_name_t {Ram_Lak_filter, Shepp_Logan_filter, Cosine_filter};
enum filter_window_t {No_window, Hann_window, Hamming_window,
		      Blackman_window};
enum filter_norm_t {CC_norm, CA_norm, PC_norm, PA_norm, SC_norm, SA_norm};

#endif // CCPI_FBP_FILTERS
