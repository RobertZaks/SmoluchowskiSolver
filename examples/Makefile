CXX = g++ -ansi -pedantic
CXXFLAGS = -Wall -I../ -O2 # -DCALC_DIRECT #-g
LIBS = -lm

BINARIES = smolsolver smoltest checktime smolsolver_obs \
	   smolgel
STATIC = 

.PHONY: default all clean

default: $(BINARIES)

all: $(BINARIES)
	cd ../smoluchowski ; $(MAKE) all

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

LIBS = -lsmoluchowski -L../smoluchowski
LIBDEPS = ../smoluchowski/libsmoluchowski.a

../smoluchowski/libsmoluchowski.a:
	cd ../smoluchowski ; $(MAKE)

smolsolver: smolsolver.cpp $(LIBDEPS)
	$(CXX) $(STATIC) $(CXXFLAGS) $^ $(LIBS) -o $@

smolsolver_obs: smolsolver_obs.cpp $(LIBDEPS)
	$(CXX) $(STATIC) $(CXXFLAGS) $^ $(LIBS) -o $@

checktime: checktime.cpp $(LIBDEPS)
	$(CXX) $(STATIC) $(CXXFLAGS) $^ $(LIBS) -o $@

smoltest: smoltest.cpp $(LIBDEPS)
	$(CXX) $(STATIC) $(CXXFLAGS) $^ $(LIBS) -o $@

smolgel: smolgel.cpp $(LIBDEPS)
	$(CXX) $(STATIC) $(CXXFLAGS) $^ $(LIBS) -o $@

clean:
	rm -rf *.o *~ $(BINARIES) *.csv times_fast times_direct
	cd scripts ; $(MAKE) clean
