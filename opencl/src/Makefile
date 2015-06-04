
CXX = g++
OCLINC = -DUSE_OPENCL -I/usr/local/cuda/include
OCLLIB = -L/usr/local/cuda/lib64 -lOpenCL
MKLINC = -I/opt/intel/composerxe/mkl/include -DMKL_ILP64
TIFFINC = -DHAS_TIFF
TIFFLIB = -ltiff
NXINC = -I/usr/local/nexus/include/nexus -DHAS_NEXUS
NXLIB = -L/usr/local/nexus/lib -lNeXusCPP
CXXFLAGS = -g -O2 -ansi -pedantic -Wall -Wno-long-long -DUSE_TIMER=true -fopenmp $(NXINC) $(TIFFINC) $(OCLINC)
LD = $(CXX)
LDFLAGS = $(CXXFLAGS)
MKLLIB = /opt/intel/composerxe/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_sequential.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_core.a
LIBS = $(NXLIB) $(TIFFLIB) $(OCLLIB)

OBJECTS = recon.o results.o instruments.o parallel.o cone.o diamond.o xtek.o\
	cgls.o p2D.o tv_reg.o tv_core.o nexus.o tiff.o timer.o utils.o \
	ui_calls.o mpi.o point.o c2D.o voxels.o landweber.o mlem.o sirt.o \
	tikhonov.o accel.o

recon: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

clean:
	rm recon *.o

$(OBJECTS): base_types.hpp allocator.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $(@F:.o=.cpp)

recon.o: recon.cpp algorithms.hpp instruments.hpp results.hpp utils.hpp \
	mpi.hpp voxels.hpp cgls.hpp tv_reg.hpp landweber.hpp mlem.hpp sirt.hpp

cgls.o: cgls.cpp algorithms.hpp instruments.hpp timer.hpp ui_calls.hpp\
	cgls.hpp blas.hpp regularize.hpp

instruments.o: instruments.cpp instruments.hpp ui_calls.hpp

parallel.o: parallel.cpp instruments.hpp parallel_f.hpp parallel_b.hpp \
	project_line.hpp timer.hpp

cone.o: cone.cpp instruments.hpp cone_f.hpp project_line.hpp timer.hpp

diamond.o: diamond.cpp instruments.hpp nexus.hpp utils.hpp ui_calls.hpp

xtek.o: xtek.cpp instruments.hpp tiff.hpp utils.hpp ui_calls.hpp

results.o: results.cpp results.hpp tiff.hpp mpi.hpp ui_calls.hpp

tv_reg.o: tv_reg.cpp algorithms.hpp instruments.hpp ui_calls.hpp blas.hpp \
	tv_reg.hpp cgls.hpp

tv_core.o: tv_core.cpp instruments.hpp algorithms.hpp blas.hpp ui_calls.hpp \
	tv_reg.hpp

landweber.o: landweber.cpp algorithms.hpp instruments.hpp timer.hpp \
	ui_calls.hpp landweber.hpp blas.hpp

mlem.o: mlem.cpp algorithms.hpp instruments.hpp timer.hpp \
	ui_calls.hpp mlem.hpp blas.hpp

sirt.o: sirt.cpp algorithms.hpp instruments.hpp timer.hpp \
	ui_calls.hpp sirt.hpp blas.hpp

tikhonov.o: tikhonov.cpp regularize.hpp

p2D.o: p2D.cpp instruments.hpp ui_calls.hpp accel.hpp

c2D.o: c2D.cpp instruments.hpp ui_calls.hpp accel.hpp

tiff.o: tiff.cpp tiff.hpp timer.hpp ui_calls.hpp

nexus.o: nexus.cpp nexus.hpp timer.hpp ui_calls.hpp

timer.o: timer.cpp timer.hpp ui_calls.hpp

utils.o: utils.cpp utils.hpp

ui_calls.o: ui_calls.cpp ui_calls.hpp

mpi.o: mpi.cpp mpi.hpp accel.hpp

accel.o: accel.cpp accel.hpp ui_calls.hpp

voxels.o: voxels.cpp voxels.hpp instruments.hpp algorithms.hpp mpi.hpp \
	utils.hpp ui_calls.hpp results.hpp blas.hpp

point.o: point.cpp
