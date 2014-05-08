
CXX = g++
MKLINC = -I/opt/intel/composerxe/mkl/include -DMKL_ILP64
TIFFINC =
TIFFLIB = -ltiff
NXINC = -I/usr/local/nexus/include/nexus -DHAS_NEXUS
NXLIB = -L/usr/local/nexus/lib -lNeXusCPP
CXXFLAGS = -g -O2 -ansi -pedantic -Wall -Wno-long-long -DUSE_TIMER=true -fopenmp $(NXINC) $(MKLINC)
LD = $(CXX)
LDFLAGS = $(CXXFLAGS)
MKLLIB = /opt/intel/composerxe/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_sequential.a /opt/intel/composerxe/mkl/lib/intel64/libmkl_core.a
LIBS = $(MKLLIB) $(NXLIB) $(TIFFLIB) 

OBJECTS = recon.o results.o instruments.o parallel.o cone.o diamond.o xtek.o\
	cgls.o tv_reg.o tv_core.o matrix2D.o nexus.o tiff.o \
	timer.o utils.o ui_calls.o mpi.o

recon: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

clean:
	rm recon *.o

$(OBJECTS): base_types.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $(@F:.o=.cpp)

recon.o: recon.cpp algorithms.hpp instruments.hpp results.hpp utils.hpp mpi.hpp

cgls.o: cgls.cpp algorithms.hpp instruments.hpp timer.hpp ui_calls.hpp

instruments.o: instruments.cpp instruments.hpp ui_calls.hpp

parallel.o: parallel.cpp instruments.hpp parallel_f.hpp parallel_b.hpp \
	project_line.hpp timer.hpp

cone.o: cone.cpp instruments.hpp cone_b.hpp cone_f.hpp \
	project_line.hpp timer.hpp ui_calls.hpp

diamond.o: diamond.cpp instruments.hpp nexus.hpp utils.hpp

xtek.o: xtek.cpp instruments.hpp tiff.hpp utils.hpp ui_calls.hpp

results.o: results.cpp results.hpp tiff.hpp mpi.hpp ui_calls.hpp

tv_reg.o: tv_reg.cpp algorithms.hpp instruments.hpp ui_calls.hpp

tv_core.o: tv_core.cpp instruments.hpp algorithms.hpp blas.hpp ui_calls.hpp

matrix2D.o: matrix2D.cpp instruments.hpp ui_calls.hpp

tiff.o: tiff.cpp tiff.hpp timer.hpp ui_calls.hpp

nexus.o: nexus.cpp nexus.hpp timer.hpp ui_calls.hpp

timer.o: timer.cpp timer.hpp ui_calls.hpp

utils.o: utils.cpp utils.hpp

ui_calls.o: ui_calls.cpp ui_calls.hpp

mpi.o: mpi.cpp mpi.hpp
