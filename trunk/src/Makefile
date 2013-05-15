
CXX = g++
CXXFLAGS = -g -O2 -ansi -pedantic -Wall -DUSE_TIMER=true -fopenmp
TIFFINC =
TIFFLIB = -ltiff
LD = $(CXX)
LDFLAGS = $(CXXFLAGS)
LIBS = $(TIFFLIB)

OBJECTS = recon.o results.o instruments.o diamond.o xtek.o cgls.o tiff.o \
	tv_reg.o tv_core.o timer.o utils.o

recon: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

clean:
	rm recon *.o

$(OBJECTS): base_types.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $(@F:.o=.cpp)

recon.o: recon.cpp algorithms.hpp instruments.hpp results.hpp utils.hpp

cgls.o: cgls.cpp algorithms.hpp instruments.hpp timer.hpp

instruments.o: instruments.cpp instruments.hpp cone_b.hpp cone_f.hpp \
	parallel_f.hpp parallel_b.hpp project_line.hpp timer.hpp

diamond.o: diamond.cpp instruments.hpp

xtek.o: xtek.cpp instruments.hpp tiff.hpp utils.hpp

results.o: results.cpp results.hpp tiff.hpp

tiff.o: tiff.cpp tiff.hpp

tv_reg.o: tv_reg.cpp algorithms.hpp instruments.hpp

tv_core.o: tv_core.cpp instruments.hpp algorithms.hpp blas.hpp

timer.o: timer.cpp timer.hpp

utils.o: utils.cpp utils.hpp
