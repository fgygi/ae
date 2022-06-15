CXX=g++
LD=$(CXX)

CXXFLAGS= $(INCLUDE) $(DFLAGS)

#LAPACKDIR=$(HOME)/software/lapack/lapack-3.5.0
#BLASIDIR=$(HOME)/software/blas/BLAS
#LIBPATH = -L$(LAPACKDIR) -L$(BLASDIR)
  
LIBS = -llapack -lblas 
 
LDFLAGS = $(LIBPATH) $(LIBS) 

all: vaepot aecheck aecheck_erf

aecheck: aecheck.o nradsolve.o vae.o
	$(LD) -o $@ $^ $(LDFLAGS)

aecheck_erf: aecheck_erf.o nradsolve.o
	$(LD) -o $@ $^ $(LDFLAGS)

vaepot: vaepot.o vae.o
	$(LD) -o $@ $^ $(LDFLAGS)

clean:
	rm -f vaepot aecheck aecheck_erf *.o
