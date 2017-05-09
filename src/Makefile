CXX=g++
LD=$(CXX)

CXXFLAGS= $(INCLUDE) $(DFLAGS)

#LAPACKDIR=$(HOME)/software/lapack/lapack-3.5.0
#BLASIDIR=$(HOME)/software/blas/BLAS
#LIBPATH = -L$(LAPACKDIR) -L$(BLASDIR)
  
LIBS = -llapack -lblas 
 
LDFLAGS = $(LIBPATH) $(LIBS) 

aecheck: aecheck.o nradsolve.o vae.o
	$(LD) -o $@ $^ $(LDFLAGS)

vaepot: vaepot.o vae.o
	$(LD) -o $@ $^ $(LDFLAGS)

clean:
	rm -f vaepot aecheck *.o 
# generate dependencies in makefile: use -Y to avoid library header files
 depend :
	makedepend -Y *.[cCh]
# DO NOT DELETE
