
CXX = g++
CXXFLAGS = -g -O2 -ansi -pedantic -Wall -fopenmp
TIFFINC =
TIFFLIB = -ltiff
LD = $(CXX)
LDFLAGS = $(CXXFLAGS)
LIBS = $(TIFFLIB)

OBJECTS = recon.o results.o instruments.o diamond.o xtek.o cgls.o tiff.o utils.o

recon: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

clean:
	rm recon *.o

$(OBJECTS): base_types.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $(@F:.o=.cpp)

recon.o: recon.cpp algorithms.hpp instruments.hpp results.hpp utils.hpp

cgls.o: cgls.cpp algorithms.hpp instruments.hpp

instruments.o: instruments.cpp instruments.hpp

diamond.o: diamond.cpp instruments.hpp parallel_b.hpp parallel_f.hpp \
	project_line.hpp

xtek.o: xtek.cpp instruments.hpp cone_b.hpp cone_f.hpp project_line.hpp \
	tiff.hpp utils.hpp

results.o: results.cpp results.hpp tiff.hpp

tiff.o: tiff.cpp tiff.hpp

utils.o: utils.cpp utils.hpp
