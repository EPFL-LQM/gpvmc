DEBUG=no
PROFILE=no
OPT=-O0 -g3
#OPT=-O2 -ffast-math -mtune=generic
#OPT=-O3 -g3 -ffast-math -funroll-loops
DEB=

INSTALL_DIR=/home/bastien/vmc-run

USEMPI=yes
STATIC=no
USEPARA=no
USEEXCEPT=no

#compilers:
INTEL=no
GNU=yes
CRAY=no
CRAY_PLAT=no

BLAS_MKL=no
BLAS_ACML=no
BLAS_ATLAS=yes
BLAS_LIBSCI=no

#random number generator
RNG_GSL=yes
RNG_MKL=no
RNG_STD=no
RNG_ACML=no

######################################################################

#computer-specifics (library paths etc...)
include local.mak

ifeq ($(RNG_GSL),yes)
    CFLAGS:=$(CFLAGS) -DUSE_RNG_GSL
else ifeq ($(RNG_MKL),yes)
    CFLAGS:=$(CFLAGS) -DUSE_RNG_MKL
else ifeq ($(RNG_STD),yes)
    CFLAGS:=$(CFLAGS) -DUSE_RNG_STD
else ifeq ($(RNG_ACML),yes)
    CFLAGS:=$(CFLAGS) -DUSE_RNG_ACML
endif

ifeq ($(INTEL),yes)
    CFLAGS:=$(CFLAGS) -DMKL_Complex8="std::complex<float>" -DMKL_Complex16="std::complex<double>"
endif

ifeq ($(PROFILE),yes)
    CFLAGS:=$(CFLAGS) -DPROFILE
endif
ifeq ($(USEPARA),yes)
    CFLAGS:=$(CFLAGS) -DUSEPARA
endif
ifeq ($(USEEXCEPT),yes)
    CFLAGS:=$(CFLAGS) -DEXCEPT
endif
ifeq ($(USEMPI),yes)
    CFLAGS:=$(CFLAGS) -DUSEMPI
endif
ifeq ($(DEBUG),yes)
    CFLAGS:=$(CFLAGS) -DDEBUG
endif

ifeq ($(STATIC),yes)
    CFLAGS:=$(CFLAGS) -static -static-libstdc++ -pthread -static-libgfortran
    LDFLAGS:=$(LDFLAGS) -static -static-libstdc++ -pthread -static-libgfortran
endif

ifeq ($(INTEL),yes)
    CFLAGS:=$(CFLAGS) -std=c++0x -DINTEL -fp-model fast=2 -U__GXX_EXPERIMENTAL_CXX0X__
    LDFLAGS:=$(LDFLAGS) -fp-model fast=2
    ifeq ($(MKL),yes)
    	CFLAGS:=$(CFLAGS) -mkl
	LDFLAGS:=$(LDFLAGS) -mkl
    endif
    ifeq ($(USEPARA),yes)
    	CFLAGS:=$(CFLAGS) -openmp -parallel
	LDFLAGS:=$(LDFLAGS) -openmp -parallel
    endif
    export OMPI_CXX=icpc
    export CXX=icpc
    export CC=icc
    export OMPI_CC=icc
    ifeq ($(USEMPI),yes)
	OCXX=mpicxx
	OFORT=mpif90
    else
    	OCXX=icpc
	OFORT=gfortran
    endif
    CFLAGS:=$(CFLAGS) $(OPT) $(DEB)
endif
ifeq ($(CRAY),yes)
    ifeq ($(GNU),no)
	CFLAGS:=$(CFLAGS) -h gnu
	CFLAGS:=$(CFLAGS) $(OPT) $(DEB)
    endif
    OCXX=CC
    OFORT=ftn
endif
ifeq ($(GNU),yes)
    CFLAGS:=$(CFLAGS) -std=c++0x -Wall -pedantic $(OPT) $(DEB)
    ifeq ($(USEPARA),yes)
    	CFLAGS:=$(CFLAGS) -fopenmp
	LDFLAGS:=$(LDFLAGS) -fopenmp
    endif
    ifeq ($(CRAY_PLAT),no)
	export OMPI_CXX=g++
	export CXX=g++
	export CC=gcc
	export OMPI_CC=gcc
	ifeq ($(USEMPI),yes)
	    OCXX=mpicxx
	    OFORT=mpif90
	else
	    OCXX=g++
	    OFORT=gfortran
	endif
    else
    	OCXX=CC
	OFORT=ftn
    endif
endif

