cmake_minimum_required(VERSION 2.8)

project(vmc CXX Fortran)

option (BuildTests "Build the tests and enable \"make test\" command." OFF)


set(VMC_SOURCES Amplitude.h SlaterDeterminant.h SlaterDeterminant.cpp ArgParse.h ArgParse.cpp BigComplex.h BigDouble.h blas_lapack.h defs.h defs.cpp FileManager.h FileManager.cpp Lattice.h LatticeState.h LatticeState.cpp LatticeStepper.h LatticeStepper.cpp linalg.h linalg.cpp MatrixQuantity.h MatrixQuantity.cpp MetroMC.h MetroMC.cpp OverlapTrack.h OverlapTrack.cpp ProjHeis.h ProjHeis.cpp Quantity.h Quantity.cpp RanGen.h RanGen.cpp ScalarQuantity.h ScalarQuantity.cpp SquareLattice.h SquareLattice.cpp StagFluxGroundState.h StagFluxGroundState.cpp StagFluxLongExciton.h StagFluxLongExciton.cpp StagFluxTransExciton.h StagFluxTransExciton.cpp StagFluxWaveFunction.h StagFluxWaveFunction.cpp StagMagnZ.h StagMagnZ.cpp StagMagnTrack.h StagMagnTrack.cpp State.h State.cpp StatSpinStruct.h StatSpinStruct.cpp Stepper.h Timer.h Timer.cpp unix_signals.h VectorQuantity.h VectorQuantity.cpp WaveFunction.h WaveFunction.cpp zdotu_sub.f)

add_executable(vmc vmc.cpp ${VMC_SOURCES})
set(VMC_TARGETS vmc)

set(VMC_INSTALL_TARGETS vmc)

option (MtROSA "Special flag to build on Monte Rosa cluster" OFF)
if (MtROSA)
    message("Build on Monte Rosa: delegate compiler flags and link to the cluster build system")
    set(CMAKE_CXX_COMPILER "CC")
    set(CMAKE_Fortran_COMPILER "ftn")
    add_definitions(-DUSEMPI)
    add_definitions(-DUSE_RNG_GSL)
else (MtROSA)
    option (USE_MPI "Use mpi to parallelize the Monte Carlo sampling" ON)
    option (USE_BLAS_ATLAS "Use ATLAS BLAS" ON)

    if (USE_BLAS_ATLAS)
        set(ENV{BLA_VENDOR} "ATLAS")
    endif (USE_BLAS_ATLAS)


    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
    find_package(HDF5 REQUIRED C HL)
    find_package(BLAS REQUIRED)
    find_package(LAPACK REQUIRED)
    find_package(FFTW REQUIRED)
    find_package(GSL)
    find_package(MPI)

    if (HDF5_FOUND)
        include_directories(${HDF5_INCLUDE_DIRS})
        foreach(t ${VMC_TARGETS})
            target_link_libraries(${t} ${HDF5_LIBRARIES})
            target_link_libraries(${t} ${HDF5_HL_LIBRARIES})
        endforeach(t)
    endif (HDF5_FOUND)
    if (BLAS_FOUND)
        include_directories(${BLAS_INCLUDE_DIRS})
        foreach(t ${VMC_TARGETS})
            target_link_libraries(${t} ${BLAS_LIBRARIES})
        endforeach(t)
    endif (BLAS_FOUND)
    if (LAPACK_FOUND)
        include_directories(${LAPACK_INCLUDE_DIRS})
        foreach(t ${VMC_TARGETS})
            target_link_libraries(${t} ${LAPACK_LIBRARIES})
        endforeach(t)
    endif (LAPACK_FOUND)
    if (FFTW_FOUND)
        include_directories(${FFTW_INCLUDE_DIRS})
        foreach(t ${VMC_TARGETS})
            target_link_libraries(${t} ${FFTW_LIBRARIES})
        endforeach(t)
    endif (FFTW_FOUND)
    if (GSL_FOUND)
        include_directories(${GSL_INCLUDE_DIRS})
        foreach(t ${VMC_TARGETS})
            target_link_libraries(${t} ${GSL_LIBRARIES})
        endforeach(t)
    endif (GSL_FOUND)

    option ( USE_RNG_GSL "Use random generator from gsl" ON )
    option ( USE_RNG_STD "Use random generator from stdlib" OFF )

    if ( USE_RNG_GSL )
        add_definitions ( -DUSE_RNG_GSL )
    endif ( USE_RNG_GSL )
    if ( USE_RNG_STD )
        add_definitions ( -DUSE_RNG_STD )
    endif ( USE_RNG_STD )

    if (USE_MPI)
        if (MPI_FOUND)
            set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
            set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
            include_directories(${MPI_INCLUDE_PATH})
            foreach(t ${VMC_TARGETS})
                target_link_libraries(${t} ${MPI_LIBRARIES})
            endforeach(t)
            add_definitions( -DUSEMPI )
        endif (MPI_FOUND)
    endif (USE_MPI)
endif (MtROSA)

option ( USE_CPP_EXCEPTION "Use C++ exception" OFF )
if (USE_CPP_EXCEPTION)
    add_definitions( -DEXCEPT )
    set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} -fexception)
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

set( CMAKE_CXX_FLAGS "-std=c++0x")

if( CMAKE_BUILD_TYPE STREQUAL "Release" )
    add_definitions( -DNDEBUG )
endif()

if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
    add_definitions( -DDEBUG )
endif()


add_subdirectory(tools)

foreach(t ${VMC_INSTALL_TARGETS})
    install(TARGETS ${t} DESTINATION bin)
endforeach(t)

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
