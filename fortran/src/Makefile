COMP=intel

ifeq ($(COMP), gnu)
    FC=mpif90  -fallow-argument-mismatch
    BLAS_ROOT=$(OPENBLAS)/lib
    BLAS=-L$(BLAS_ROOT) -lopenblas
    LAPACK=
    OMP=-fopenmp
endif
ifeq ($(COMP), intel)
    FC=mpiifort
    MKLPATH=$(MKLROOT)/lib/intel64
    MKLINCLUDE=$(MKLROOT)/include
    BLAS=-L$(MKLPATH) -I$(MKLINCLUDE) -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm
    LAPACK=-L$(MKLPATH) -I$(MKLINCLUDE) -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm
    OMP=-qopenmp
endif
ifeq ($(COMP), nvidia)
    FC=mpif90
    LAPACK=-llapack
    BLAS=-cudalib=cublas -lblas
    OMP=-mp
endif
ifeq ($(COMP), amd)
    FC=mpif90
    AMDLIB_ROOT=$(AOCLLIB)/lib
    AMDBLIS=$(AMDLIB_ROOT)
    AMDLIBFLAME=$(AMDLIB_ROOT)
    BLAS=-L$(AMDBLIS) -lblis-mt
    LAPACK=-L$(AMDLIBFLAME) -lflame
    OMP=-fopenmp
endif

LIBS=$(BLAS) $(LAPACK)

default: linpack_test.o

linpack.o: linpack.f90
	$(FC) $(OMP) -c linpack.f90
linpack_test.o: linpack_test.f90 linpack.o
	$(FC) $(OMP) -o linpack_test linpack_test.f90 *.o $(LIBS)
clean:
	rm -f *.o *.mod linpack_test
