Variational Monte-Carlo for Gutzwiller-Projected Wavefunctions implementation {#mainpage}
=============================================================================

This software was used to produce the numerical results published in

B. Dalla Piazza _et al._, "Fractional excitations in the square lattice quantum antiferromagnet", Accepted in Nature Physics (October 2014).

The software main objective is to provide an efficient way of calculating many-body wavefunction $|\psi\rangle$ overlaps with real space spin configuration state $|\alpha\rangle$. Its main components are:

*   The wavefunction $|\psi\rangle$ is a trial ground state or excited state. Different such wavefunctions are implemented and inherit from the WaveFunction object.
*   The real space spin configuration $|\alpha\rangle$ are states upon which a Monte Carlo sampling is performed. They are implemented in the LatticeState object.
*   The overlap $\langle\psi|\alpha\rangle$ is a Slater determinant and is efficiently implemented into the SlaterDeterminant object.
*   The Monte Carlo random walk generator (_i. e._ the step generator and the Monte Carlo weight definition) is implemented in the objects inheriting from the Stepper object.
*   Different measurements can be performed simultaneously during the Monte Carlo sampling. These inherit from the Quantity object.
*   Finally the Monte Carlo itself is implemented into the MetroMC object (for Metropolis Monte Carlo).

##Requirements
This software has been built in various linux distributions. The requirements are:
*   A BLAS/LAPACK library
*   The HDF5 library
*   cmake

##Optional requirements
*   The GNU Scientific Library
*   An MPI implementation
*   doxygen to build the documentation

##Install

1.  Make a build directory, for instance directly in the main source folder
        :~/.../sourcedir$ mkdir build
2.  Configure the build toolchain.
        :~/.../sourcedir$ cd build
        :~/.../sourcedir/build$ cmake ..
3.  Review the options to check they meet your environment
        :~/.../sourcedir/build$ ccmake .
In particular:
    *   If you don't have the GNU Scientific Library, please choose "USE_RNG_GSL=NO" and "USE_RNG_STD=YES". Doing so you will use the standard random number generator (probably a bad idea).
    *   If you don't have an MPI implementation installed, you must pick "USE_MPI=NO".
    *   Set the installation directory (default /usr/local) to some place you have write access (or you can do sudo)
4.  Build the software and install it
        :~/.../sourcedir/build$ make
        :~/.../sourcedir/build$ make install
This produces the "vmc" executable and installs python packages (in the installation directory) that you can use for postprocessing (Make sure python knows where to find them).

##Run
The "vmc" executable, if called without arguments, will list all the possible options. A typical example would be:

        :~/.../installdir$ vmc --neel=0.1 --phi=0.25 --L=8\
        --samples=1 --samples_saves=1 --samples_saves_stat=100000 --therm=200\
        --channel=groundstate --meas_projheis --meas_stagmagn --meas_statspinstruct

*   --neel and --phi options relate to some variational parameters (see our publication above)
*   --L=8 sets the system size to 8x8
*   --samples=1 sets the number of different samples (per thread) to gather statistics for.
*   --samples\_saves=1 sets the number of times the gathered statistics must be saved.
*   --samples\_saves\_stat=100000 sets the amount of statistics of each save.
*   --therm=200 sets the thermalization steps to 200\*L\*L
*   --channel=groundstate asks to measure the groundstate trial wavefunction.
*   --meas\_projheis asks to measure the (projected) Heisenberg Hamiltonian.
*   --meas\_stagmagn asks to measure the staggered magnetization.
*   --meas\_statspinstruct asks to measure the static spin structure factor.

A python script performing an example calculation may be found in the source root directory.

##Output
The different outputs are in the HDF5 format, one file per quantity. The HDF5 files gather all the statistics of each quantity. Some postprocessing tools are provided in python.

##License
The Software is distributed under the MIT license.
