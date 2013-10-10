################################################################################
#
# Build options for xray project
#
# Optional command-line flags to control build options:
#
# dbg       - Set to 1 to build a debug version. Default: unset.
# hmxif     - Set to 1 to build code with _HMXIF defined. Default: unset
# verbose   - Set to 1 to echo compiler calls. Default unset.
# fastmath  - Set to 1 to add -use_fast_math to nvcc (not recommended)
# nvcc_warn_verbose - Set to 1 to pass through CXXWARN_FLAGS to nvcc
# maxregisters      - Set to the max number of registers to allow cuda code to use
#
################################################################################

# Set paths specific to your system (here or in your environment)
CUDA_INSTALL_PATH     ?= /opt/cuda/3.2/cuda
CUDA_SDK_INSTALL_PATH ?= $(HOME)/dev/nvidia/cuda/3.2/NVIDIA_GPU_Computing_SDK
INTEL_INSTALL_PATH    ?= /opt/intel/Compiler/11.0/074

# Add source files here
EXECUTABLE	:= cuFBP
# Cuda source files (compiled with cudacc)
CUFILES		:= \
 cone_gpu.cu \
 parallel_gpu.cu
# CUDA dependency files
CU_DEPS		:= defs_gpu.cuh
# C/C++ source files
# Upper half of this list is the set of files required
# for a Diamond build. The lower half of additional files
# are needed for an HMXIF build. We can safely include
# all files for all builds due to conditional source code
# compilation via the _HMXIF macro in source the files.
CCFILES	:= \
 allFDK.cpp \
 bp_all.cpp \
 diamond_cpp.cpp \
 flatDark.cpp \
 gen_functions.cpp \
 gtiff.cpp \
 HMutil.cpp \
 main.cpp \
 MatPar.cpp \
 message.cpp \
 newData.cpp \
 ole2cd.cpp \
 parallel.cpp \
 PolarToCart.cpp \
 readWriteSlice.cpp \
 ringArtefacts.cpp \
 time_stamp.cpp \
 transformSinogram.cpp \
 xmlStructure.cpp \
 xrm2tif.cpp \
 xtiff.cpp

# Additional libraries needed by the project
USECUFFT        := 1

################################################################################
# Rules and targets

include common.mk
