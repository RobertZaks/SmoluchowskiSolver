CXX = g++ -ansi -pedantic
CXXFLAGS = -Wall -O2 #-g
CPPFLAGS = -D HANDLEXCEPTION=0
LIBS = -lm

FILES = matrix.o indexlist.o skeleton.o cross.o fft.o operator.o solver.o

LIBNAME = smoluchowski
LIBFILES = lib$(LIBNAME).a


lib$(LIBNAME).a: $(FILES)
	$(AR) -crs $@ $(FILES)

.PHONY: all clean

all: lib$(LIBNAME).a testcross testindexlist testfft

%.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

testfft: testfft.o fft.o
	$(CXX) $(STATIC) $(CXXFLAGS) $(CPPFLAGS) $^ $(LIBS) -o $@

testindexlist: testindexlist.o indexlist.o
	$(CXX) $(STATIC) $(CXXFLAGS) $^ $(LIBS) -o $@

testcross: testcross.o matrix.o indexlist.o skeleton.o cross.o
	$(CXX) $(STATIC) $(CXXFLAGS) $(CPPFLAGS) -D CALCNORMS=1 $^ $(LIBS) -o $@

clean:
	rm -f *.o *~ testfft testcross testindexlist lib$(LIBNAME).a
