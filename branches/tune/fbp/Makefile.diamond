#!/usr/local/bin/gnumake


#CXX := icpc
CXX := g++
LINK := g++
#CXX := g++
#CXX := icpc
.SUFFIXES : .cpp .cu .o 

maxregisters := 16




#INTELPREFIX := /home/kny48981/local/intel/

INTELPREFIX := /dls_sw/apps/intel/
CUDAPREFIX := /usr/local/cuda/
CUDASDKPREFIX := /dls_sw/apps/NVIDIA_GPU_Computing_SDK/C
NC := $(CUDAPREFIX)/bin/nvcc

OPTIMIZER = -O3
#OPTIMIZER = -g -DGNU_TRAP_FPE
#OPTIMIZER = -g 

#exclude all warnings ..  Too many warnings!
#FLAGS = -w $(OPTIMIZER) 

#include all warnings .. 
FLAGS = -Wall $(OPTIMIZER) 

FLAGS += -DDIAMOND # -D_ALGORITHM1
#FLAGS = -w -ggdb
FLAGS += $(CFLAGS)

ifdef maxregisters
	NVCCFLAGS += -maxrregcount=$(maxregisters)
endif
NVCCFLAGS += -use_fast_math

#location of the library folders

LIBDIRS = -L$(CUDASDKPREFIX)/lib  \
-L$(CUDAPREFIX)/lib64/ \
-L$(INTELPREFIX)/ipp/em64t/sharedlib \
-L${INTELPREFIX}/mkl/lib/em64t \
-L${INTELPREFIX}/tbb/em64t/cc4.1.0_libc2.4_kernel2.6.16.21/lib \
-L${INTELPREFIX}/lib/intel64/ \
-L/lib64 

INCDIRS := \
-I$(CUDAPREFIX)/include/ \
-I${INTELPREFIX}/ipp/include \
-I$(CUDASDKPREFIX)/common/inc/

#library flags required
#
#


LIBS = -lm   -ldl -lrt \
-lcutil \
-lcudart\
-lcufft \
-Wl,-rpath $(CUDAPREFIX)/lib \
-Wl,-rpath $(CUDAPREFIX)/lib64 \
-Wl,-rpath ${INTELPREFIX}/ipp/em64t/sharedlib \
-Wl,-rpath ${INTELPREFIX}/mkl/lib/em64t \
-Wl,-rpath ${INTELPREFIX}/tbb/em64t/cc4.1.0_libc2.4_kernel2.6.16.21/lib \
-Wl,-rpath ${INTELPREFIX}/lib/intel64/ \
-Wl,-rpath-link ${INTELPREFIX}/ipp/em64t/sharedlib \
-Wl,-rpath-link ${INTELPREFIX}/mkl/lib/em64t \
-Wl,-rpath-link ${INTELPREFIX}/tbb/em64t/cc4.1.0_libc2.4_kernel2.6.16.21/lib \
-Wl,-rpath-link ${INTELPREFIX}/lib/intel64/ \
-Wl,-rpath ${INTELPREFIX}lib  \
-Wl,-rpath $(CUDASDKPREFIX)/lib  \
-pthread \
-lippcoreem64t \
-lippcoreem64t \
-lippcvem64t \
-lippcvem64t \
-lippcvu8 \
-lippchem64t \
-lippchem64t \
-lippchu8 \
-lippiem64t \
-lippiem64t \
-lippiu8 \
-lippsem64t \
-lippsem64t \
-lippsu8 \
-lippvmem64t \
-lippvmem64t \
-lippvmu8 \
-lmkl_core \
-lmkl_intel_lp64 \
-lmkl_intel_thread \
-lmkl_vml_def \
-lmkl_vml_mc \
-ltiff \
-ltbb 

#
#-lgomp\
#
#
#-Wl,--verbose \

CUDASOURCES=\
parallel_gpu.cu \
cone_gpu.cu

CSOURCES=\
allFDK.cpp \
bp_all.cpp \
diamond_cpp.cpp\
flatDark.cpp \
gen_functions.cpp \
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
gtiff.cpp\
handlers.cpp\
xtiff.cpp 


# I think you could use a smaller set of source files
CSOURCES_MINIMAL=\
bp_all.cpp \
diamond_cpp.cpp \
gtiff.cpp \
xtiff.cpp \
newData.cpp \
xmlStructure.cpp \
message.cpp \
gen_functions.cpp \
main.cpp \
parallel.cpp \
readWriteSlice.cpp \
ringArtefacts.cpp \
flatDark.cpp \
transformSinogram.cpp \
xrm2tif.cpp \
time_stamp.cpp

COBJECTS=$(CSOURCES:%.cpp=%.o)
CUDAOBJECTS=$(CUDASOURCES:%.cu=%.o)

PROGNAME=dt64n-alg1




# the dependencies
default: $(PROGNAME)

.cpp.o:
	$(CXX) $(INCDIRS) $(FLAGS) -c $<

.cu.o:
	$(NC) $(INCDIRS) $(OPTIMIZER) -c $< $(NVCCFLAGS)


$(PROGNAME): $(CUDAOBJECTS) $(COBJECTS) 
	$(CXX) -o $@ $(CUDAOBJECTS) $(COBJECTS) $(LIBDIRS) $(LIBS) $(FLAGS)

# using 'export NVCCFLAGS=-DDEBUGTIME' in the shell
# to activate the extra time print out


#handlers.o:	handlers.cpp
#	$(CXX) -w $(INCDIRS) $(FLAGS) -c handlers.cpp


cclean:
	rm  -f $(CUDAOBJECTS)  $(PROGNAME)

clean:
	rm  -f $(COBJECTS) $(CUDAOBJECTS)  $(PROGNAME)


