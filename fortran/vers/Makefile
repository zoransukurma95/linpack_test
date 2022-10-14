COMP=intel

ifeq ($(COMP), gnu)
    FC=mpif90  -fallow-argument-mismatch
    BLAS_ROOT=/gpfs/opt/sw/spack-0.17.1/opt/spack/linux-almalinux8-zen3/gcc-11.2.0/netlib-lapack-3.9.1-y3qzs3hekhmfhr2g5xpzh63ijgrkivno/lib64
    LAPACK_ROOT=/gpfs/opt/sw/spack-0.17.1/opt/spack/linux-almalinux8-zen3/gcc-11.2.0/netlib-lapack-3.9.1-y3qzs3hekhmfhr2g5xpzh63ijgrkivno/lib64
    BLAS=$(BLAS_ROOT)/libblas.a
    LAPACK=$(LAPACK_ROOT)/liblapack.a
endif
ifeq ($(COMP), intel)
    FC=mpiifort
    MKLPATH=$(MKLROOT)/lib/intel64
    MKLINCLUDE=$(MKLROOT)/include
    BLAS=-L$(MKLPATH) -I$(MKLINCLUDE) -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm
    LAPACK=-L$(MKLPATH) -I$(MKLINCLUDE) -Wl,--start-group $(MKLPATH)/libmkl_intel_lp64.a $(MKLPATH)/libmkl_intel_thread.a $(MKLPATH)/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm
endif
ifeq ($(COMP), nvidia)
    FC=mpif90
    LAPACK=-llapack
    BLAS=-cudalib=cublas -lblas
endif
ifeq ($(COMP), amd)
    FC=mpif90
    AMDLIB_ROOT=$(HOME)/lib/AMD-AOCL/3.2.0/lib
    AMDBLIS=$(AMDLIB_ROOT)
    AMDLIBFLAME=$(AMDLIB_ROOT)
    BLAS=$(AMDBLIS)/libblis.a
    LAPACK=$(AMDLIBFLAME)/libflame.a
endif

LIBS=$(BLAS) $(LAPACK)

default: linpack_test.o

linpack.o: linpack.f90
	$(FC) -c linpack.f90
linpack_test.o: linpack_test.f90 linpack.o
	$(FC) -o linpack_test linpack_test.f90 *.o $(LIBS)
clean:
	rm -f *.o *.mod linpack_test
