cmake_minimum_required(VERSION 2.8)

project(vmc CXX Fortran)

add_executable(vmc Amplitude.h Amplitude.cpp ArgParse.h ArgParse.cpp BigComplex.h BigDouble.h blas_lapack.h defs.h defs.cpp FileManager.h FileManager.cpp Lattice.h LatticeState.h LatticeState.cpp LatticeStepper.h LatticeStepper.cpp linalg.h linalg.cpp Magnetization.h Magnetization.cpp MatrixQuantity.h MatrixQuantity.cpp MetroMC.h MetroMC.cpp OverlapTrack.h OverlapTrack.cpp ProjHeis.h ProjHeis.cpp Quantity.h Quantity.cpp RanGen.h RanGen.cpp ScalarQuantity.h ScalarQuantity.cpp SFpNpHxExciton.h SFpNpHxExciton.cpp SFpNpHxGroundState.h SFpNpHxGroundState.cpp SFpNpHxWaveFunction.h SFpNpHxWaveFunction.cpp SquareLattice.h SquareLattice.cpp StagFluxGroundState.h StagFluxGroundState.cpp StagFluxLongExciton.h StagFluxLongExciton.cpp StagFluxTransExciton.h StagFluxTransExciton.cpp StagFluxWaveFunction.h StagFluxWaveFunction.cpp StagMagn.h StagMagn.cpp StagMagnTrack.h StagMagnTrack.cpp State.h State.cpp StatSpinStruct.h StatSpinStruct.cpp Stepper.h Timer.h Timer.cpp unix_signals.h VectorQuantity.h VectorQuantity.cpp vmc.cpp WaveFunction.h WaveFunction.cpp zdotu_sub.f)

option (USE_MPI "Use mpi to parallelize the Monte Carlo sampling" ON)
option (USE_BLAS_ATLAS "Use ATLAS BLAS" ON)
option (USE_BLAS_SCSL "Use SGI/Cray BLAS" OFF)

if (USE_BLAS_ATLAS)
    set(ENV{BLA_VENDOR} "ATLAS")
elseif (USE_BLAS_SCSL)
    set(ENV{BLA_VENDOR "SCSL")
endif(USE_BLAS_ATLAS)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(HDF5 REQUIRED C HL)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(GSL)
find_package(MPI)

if (HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIRS})
    target_link_libraries(vmc ${HDF5_LIBRARIES})
    target_link_libraries(vmc ${HDF5_HL_LIBRARIES})
endif (HDF5_FOUND)
if (BLAS_FOUND)
    include_directories(${BLAS_INCLUDE_DIRS})
    target_link_libraries(vmc ${BLAS_LIBRARIES})
endif (BLAS_FOUND)
if (LAPACK_FOUND)
    include_directories(${LAPACK_INCLUDE_DIRS})
    target_link_libraries(vmc ${LAPACK_LIBRARIES})
endif (LAPACK_FOUND)
if (GSL_FOUND)
    include_directories(${GSL_INCLUDE_DIRS})
    target_link_libraries(vmc ${GSL_LIBRARIES})
    add_definitions( -DUSE_RNG_GSL )
else (GSL_FOUND)
    add_definitions( -DUSE_RNG_STD )
endif (GSL_FOUND)

if (USE_MPI)
    if (MPI_FOUND)
        set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
        set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
        include_directories(${MPI_INCLUDE_PATH})
        target_link_libraries(vmc ${MPI_LIBRARIES})
        add_definitions( -DUSEMPI )
    endif (MPI_FOUND)
endif (USE_MPI)

option ( USE_CPP_EXCEPTION "Use C++ exception" OFF )
if (USE_CPP_EXCEPTION)
    add_definitions( -DEXCEPT )
endif (USE_CPP_EXCEPTION)

option ( USE_PROFILING "Turn on profiling of the code (timers)" OFF )
if ( USE_PROFILING )
    add_definitions( -DPROFILE )
endif( USE_PROFILING )

exec_program(
    "git"
    ${CMAKE_CURRENT_SOURCE_DIR}
    ARGS "describe"
    OUTPUT_VARIABLE VERSION )

string( REGEX MATCH "-g.*$" VERSION_SHA1 ${VERSION} )
string( REGEX REPLACE "[-g]" "" VERSION_SHA1 ${VERSION_SHA1} )

add_definitions( -DGIT_SHA1="${VERSION_SHA1}" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -ffast-math -funroll-loops -march=native" CACHE STRING "Compile flags" FORCE)

install(TARGETS vmc DESTINATION .)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
        )
endif(DOXYGEN_FOUND)
