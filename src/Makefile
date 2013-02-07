
CXX = g++
CXXFLAGS = -g -ansi -pedantic -Wall -fopenmp
LD = $(CXX)
LDFLAGS = $(CXXFLAGS)
LIBS =

OBJECTS = recon.o results.o diamond.o xtek.o cgls.o

recon: $(OBJECTS)
	$(LD) $(LDFLAGS) -o $@ $(OBJECTS) $(LIBS)

clean:
	rm recon *.o

$(OBJECTS): base_types.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $(@F:.o=.cpp)

recon.o: recon.cpp algorithms.hpp instruments.hpp results.hpp

cgls.o: cgls.cpp algorithms.hpp instruments.hpp

diamond.o: diamond.cpp instruments.hpp

xtek.o: xtek.cpp instruments.hpp xtek_b.hpp xtek_f.hpp project_line.hpp

results.o: results.cpp results.hpp
