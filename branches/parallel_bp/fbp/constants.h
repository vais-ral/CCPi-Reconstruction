#define MAX_MESSAGE 1024

#define warp_size 32
#define shared_elements 1024

#define pi  3.14159265358979
#define pif  3.14159265358979f

#define nw 1
#define ncos 10
#define gpu_filter_chunk 16
#define block_height 12
#define block_width 16
#define line_width 1000
#define name_width 100
#define value_width 1000
#define len_name 50

#define b_s 1
#define b_r 16
#define b_a 32

#define RAW_TYPE_NO 0
#define RAW_TYPE_YES 1

#define BYTE_ORDER_BIG 0
#define BYTE_ORDER_LITTLE 1

#define SHAPE_POINT 0
#define	SHAPE_PIXEL 1
#define SHAPE_INPUT 2

#define SHIFTS_NO 0
#define SHIFTS_POSITIVE 1
#define SHIFTS_NEGATIVE 2

#define SHIFTS_ABSOLUTE 5
#define SHIFTS_RELATIVE 6

#define NAME_FBP 0
#define NAME_XRM2TIF 1
#define NAME_ALLFDK 2

#define XRM_OUTPUT_ALL 0 
#define	XRM_OUTPUT_XML 1

#define MAX_ERROR 1024
#define MAX_FOLDER 512
#define MAX_FILE 200
#define MAX_EXT  10
#define MAX_CHILD 50
#define MAX_ENV 100

#define INPUT_TYPE_INTENSITY 0
#define INPUT_TYPE_ATTENUATION 1

#define INPUT_RESTRICTIONS_NO 0
#define INPUT_RESTRICTIONS_YES 1

#define INPUT_ORIENTATION_SLICE 0
#define INPUT_ORIENTATION_PROJECTION 1

#define FIELD_PROFILE_CONSTANT 0
#define FIELD_PROFILE_LINEAR 1
#define FIELD_PROFILE_PROFILE 3

#define FIELD_TYPE_USER 0
#define FIELD_TYPE_ROW 5

#define FIELD_TYPE_VALUE 6
#define FIELD_TYPE_FILE 7

#define HIGH_PEAKS_NO 0
#define HIGH_PEAKS_YES 1

#define INTENSITY_NO 0
#define INTENSITY_ROW 1
#define INTENSITY_COLUMN 2
#define INTENSITY_ZERO 3
 
#define RING_ARTEFACTS_NO 0
#define RING_ARTEFACTS_COLUMN 1
#define RING_ARTEFACTS_AML 2

#define MISSED_PROJECTIONS_NO 0
#define MISSED_PROJECTIONS_LINEAR 1
#define MISSED_PROJECTIONS_ZERO 2

#define ROTATION_ANGLE_180 0
#define ROTATION_ANGLE_OTHER 1

#define ROTATION_ANGLE_END_NO 0
#define ROTATION_ANGLE_END_YES 1

#define SCALE_NO 0
#define SCALE_YES 1

#define EXTRAPOLATION_NO 0
#define EXTRAPOLATION_ZERO 1
#define EXTRAPOLATION_CONSTANT 2


#define CLOCKWISE_ROTATION_NO 0
#define CLOCKWISE_ROTATION_YES 1
		
#define FILTER_NO 0
#define FILTER_YES 1

#define FILTER_NAME_RL 0
#define FILTER_NAME_SL 1
#define FILTER_NAME_COS 2

#define FILTER_WINDOW_NO 0
#define FILTER_WINDOW_HANN 1
#define FILTER_WINDOW_HAMMING 2
#define FILTER_WINDOW_BLACKMAN 3

#define FILTER_NORM_CC 0
#define FILTER_NORM_CA 1
#define FILTER_NORM_PC 2
#define FILTER_NORM_PA 3
#define FILTER_NORM_SC 4
#define FILTER_NORM_SA 5

#define TILT_NO 0
#define TILT_YES 1


#define RADIUS_PIXEL 0
#define RADIUS_PERCENT 1

#define PC_INTERPOLATION_NN 0
#define PC_INTERPOLATION_LINEAR 1
#define PC_INTERPOLATION_CUBIC 2


#define ROI_STANDARD 0
#define ROI_RECTANGLE 1


#define OUTPUT_WIDTH_STANDARD 0
#define OUTPUT_WIDTH_GIVEN 1

#define OUTPUT_TYPE_SOLUTION 0
#define OUTPUT_TYPE_FFCORRECTION 1
#define OUTPUT_TYPE_PREPROCESSED 2
#define OUTPUT_TYPE_TRANSFORMED 3

#define OUTPUT_STATE_ATTENUATION 0
#define OUTPUT_STATE_INTENSITY 1

#define OUTPUT_BITS_INPUT 0
#define OUTPUT_BITS_GIVEN 1

#define OUTPUT_RESTRICTIONS_NO 0
#define OUTPUT_RESTRICTIONS_YES 1

#define READ_DARK 0
#define READ_FLAT 1

#define MAX_TAG_LENGTH 50
#define MAX_TAG_COUNT 10

#define XRM_UNKNOWN -1 
#define XRM_NOI 0 
#define XRM_HEIGHT 1
#define XRM_WIDTH 2
#define XRM_PIXEL_SIZE 3
#define XRM_DATA_TYPE 4
#define XRM_D2A 5
#define XRM_S2A 6


#define DEFAULT_GPU_BLOCK_WIDTH 0
#define DEFAULT_GPU_BLOCK_HEIGHT 0
