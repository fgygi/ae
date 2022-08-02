CXX=g++
LD=$(CXX)

CXXFLAGS= $(INCLUDE) $(DFLAGS)

#LAPACKDIR=$(HOME)/software/lapack/lapack-3.5.0
#BLASIDIR=$(HOME)/software/blas/BLAS
#LIBPATH = -L$(LAPACKDIR) -L$(BLASDIR)
FFTWDIR=$(HOME)/software/fftw/fftw-3.3.4
  
LIBS = -llapack -lblas 
 
LDFLAGS = $(LIBPATH) $(LIBS) 

all: vaepot aecheck aecheck_erf

aecheck: aecheck.o nradsolve.o vae.o
	$(LD) -o $@ $^ $(LDFLAGS)

aecheck_erf: aecheck_erf.o nradsolve.o
	$(LD) -o $@ $^ $(LDFLAGS)

vaepot: vaepot.o vae.o
	$(LD) -o $@ $^ $(LDFLAGS)

sinft.o: sinft.cpp
	$(CXX) -I$(FFTWDIR)/api -DUSE_FFTW3 -c $^

vaefour: vaefour.o sinft.o vae.o
	$(LD) -o $@ -I$(FFTWDIR)/api -DUSE_FFTW3 $^ -L$(FFTWDIR)/.libs -lfftw3

gfour: gfour.o sinft.o
	$(LD) -o $@ -I$(FFTWDIR)/api -DUSE_FFTW3 $^ -L$(FFTWDIR)/.libs -lfftw3


clean:
	rm -f vaepot aecheck aecheck_erf *.o
