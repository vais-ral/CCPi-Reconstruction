
CXX = g++
PYTHONINC = -I/usr/include/python -I/usr/lib64/python2.7/site-packages/numpy/core/include
PYTHONLIB = -lboost_python
TIFFINC =
TIFFLIB =
NXINC = -I/usr/local/nexus/include/nexus -DHAS_NEXUS
NXLIB = -L/usr/local/nexus/lib -lNeXusCPP
CXXFLAGS = -g -O2 -fPIC -ansi -Wall -Wno-long-long -DUSE_TIMER=true -fopenmp -I../src $(PYTHONINC)
LD = $(CXX)
LDFLAGS = $(CXXFLAGS) -shared
LIBS = $(PYTHONLIB)

PYOBJS = ccpi_wrapper.o
OBJECTS = results.o instruments.o parallel.o cone.o diamond.o xtek.o cgls.o \
	p2D.o tv_reg.o tv_core.o tiff.o timer.o utils.o ui_calls.o \
	mpi.o c2D.o voxels.o sirt.o mlem.o tikhonov.o total_v.o accel.o

reconstruction.so: $(PYOBJS) $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(PYOBJS) $(OBJECTS) $(LIBS)

clean:
	rm reconstruction.so *.o

$(PYOBJS): ../src/base_types.hpp ../src/allocator.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $(@F:.o=.cpp)

$(OBJECTS): ../src/base_types.hpp ../src/allocator.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ ../src/$(@F:.o=.cpp)

ccpi_wrapper.o: ccpi_wrapper.cpp reconstruct.hpp reconstruct.cpp \
	../src/mpi.hpp ../src/utils.hpp ../src/instruments.hpp \
	../src/algorithms.hpp ../src/results.hpp ../src/voxels.hpp \
	../src/cgls.hpp ../src/sirt.hpp ../src/mlem.hpp

cgls.o: ../src/cgls.cpp ../src/algorithms.hpp ../src/instruments.hpp \
	../src/timer.hpp ../src/ui_calls.hpp ../src/blas.hpp ../src/cgls.hpp \
	../src/regularize.hpp

mlem.o: ../src/mlem.cpp ../src/algorithms.hpp ../src/instruments.hpp \
	../src/timer.hpp ../src/ui_calls.hpp ../src/mlem.hpp ../src/blas.hpp

sirt.o: ../src/sirt.cpp ../src/algorithms.hpp ../src/instruments.hpp \
	../src/timer.hpp ../src/ui_calls.hpp ../src/sirt.hpp ../src/blas.hpp

tikhonov.o: ../src/tikhonov.cpp ../src/regularize.hpp

total_v.o: ../src/total_v.cpp ../src/regularize.hpp

total_v.o: ../src/total_v.cpp ../src/regularize.hpp

instruments.o: ../src/instruments.cpp ../src/instruments.hpp ../src/ui_calls.hpp

parallel.o: ../src/parallel.cpp ../src/instruments.hpp ../src/parallel_f.hpp \
	../src/parallel_b.hpp ../src/project_line.hpp ../src/timer.hpp

cone.o: ../src/cone.cpp ../src/instruments.hpp ../src/cone_f.hpp \
	../src/project_line.hpp ../src/timer.hpp

diamond.o: ../src/diamond.cpp ../src/instruments.hpp \
	../src/utils.hpp ../src/ui_calls.hpp

xtek.o: ../src/xtek.cpp ../src/instruments.hpp ../src/tiff.hpp \
	../src/utils.hpp ../src/ui_calls.hpp

results.o: ../src/results.cpp ../src/results.hpp ../src/tiff.hpp \
	../src/mpi.hpp ../src/ui_calls.hpp

tv_reg.o: ../src/tv_reg.cpp ../src/algorithms.hpp ../src/instruments.hpp \
	../src/ui_calls.hpp ../src/blas.hpp ../src/cgls.hpp ../src/tv_reg.hpp

tv_core.o: ../src/tv_core.cpp ../src/instruments.hpp ../src/algorithms.hpp \
	../src/blas.hpp ../src/ui_calls.hpp ../src/tv_reg.hpp

p2D.o: ../src/p2D.cpp ../src/instruments.hpp ../src/ui_calls.hpp

c2D.o: ../src/c2D.cpp ../src/instruments.hpp ../src/ui_calls.hpp

tiff.o: ../src/tiff.cpp ../src/tiff.hpp ../src/timer.hpp ../src/ui_calls.hpp

timer.o: ../src/timer.cpp ../src/timer.hpp ../src/ui_calls.hpp

utils.o: ../src/utils.cpp ../src/utils.hpp

ui_calls.o: ../src/ui_calls.cpp ../src/ui_calls.hpp

mpi.o: ../src/mpi.cpp ../src/mpi.hpp ../src/accel.hpp

accel.o: ../src/accel.cpp ../src/accel.hpp ../src/ui_calls.hpp

voxels.o: ../src/voxels.cpp ../src/instruments.hpp ../src/utils.hpp \
	../src/voxels.hpp ../src/algorithms.hpp ../src/mpi.hpp \
	../src/ui_calls.hpp ../src/results.hpp ../src/blas.hpp
