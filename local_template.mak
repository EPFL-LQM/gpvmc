# MKL Compile and link flags
MKLINCLUDE=-I<intel mkl include path>
MKLLINK=<look up at http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/>

# ACML Compile and link flags
ACMLINCLUDE=-I<acml include path>
ACMLLINK=-L<acml library path> -lacml<optional: _mp, _fma4, _fma4_mp>

# ATLAS compile and link flags
ATLASINCLUDE=-I/usr/include/atlas
ATLASLINK=-llapack

# GSL compile and link flags
GSLINCLUDE=-I/usr/include/gsl
GSLLINK=-lgsl

# HDF5 compile and link flags
HDF5INCLUDE=
HDF5LINK=-lhdf5 -lhdf5_hl

ifeq ($(BLAS_MKL),yes)
    CFLAGS:= $(CFLAGS) $(MKLINCLUDE)
    ifeq ($(GNU),yes)
    	CFLAGS:=$(CFLAGS) -m64
	LDFLAGS:=$(LDFLAGS) $(MKLLINK)
    endif
else ifeq ($(RNG_MKL),yes)
    CFLAGS:= $(CFLAGS) $(MKLINCLDE)
    ifeq ($(GNU),yes)
    	CFLAGS:=$(CFLAGS) -m64
	LDFLAGS:=$(LDFLAGS) $(MKLLINK)
    endif
endif

# Options for AMD ACML
ifeq ($(BLAS_ACML),yes)
    CFLAGS:=$(CFLAGS) $(ACMLINCLUDE) -m64
    LDFLAGS:=$(LDFLAGS) $(ACMLLINK)
else ifeq ($(RNG_ACML),yes)
    CFLAGS:=$(CFLAGS) $(ACMLINCLUDE) -m64
    LDFLAGS:=$(LDFLAGS) $(ACMLLINK)
endif

# Options for ATLAS blas
ifeq ($(BLAS_ATLAS),yes)
    CFLAGS:= $(CFLAGS) $(ATLASINCLUDE)
    LDFLAGS:=$(LDFLAGS) $(ATLASLINK)
endif

# Options for the GSL random number generator
ifeq ($(RNG_GSL),yes)
    CFLAGS:=$(CFLAGS) $(GSLINCLUDE)
    LDFLAGS:=$(LDFLAGS) $(GSLLINK)
endif

CFLAGS:=$(CFLAGS) $(HDF5INCLUDE)
LDFLAGS:=$(LDFLAGS) $(HDF5LINK)
