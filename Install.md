
# Vulture Finite Difference Time-Domain Electromagnetic Solver


## Quick start for Linux

    git clone <path-to-repository> vulture-working
    mkdir vulture-build-linux_$(arch)
    cd vulture-build-linux_$(arch)
    cmake -D CMAKE_INSTALL_PREFIX=/usr/local -D CMAKE_BUILD_TYPE=Release -D WITH_OPENMP=ON \
       -D USE_SCALED_FIELDS=ON -D WITH_TESTS=ON -D GNUPLOT_TESTS=ON -D PYTHON_TESTS=ON \
       -D OCTAVE_TESTS=ON -D PROCESSING_TESTS=ON -D BIG_TESTS=OFF -D WITH_SIBC=ON ../vulture-working/
    make
    make test
    make install
    make package
    make package_source

## Building from source on Linux

Requirements: 

 * [git][], optional, for obtaining source from repository
 
 * [CMake][], required, version >= 2.8.3
 
 * C compiler, , required, with C99 support
 
 * [Octave][] or [MATLAB][], optional, for some more advanced tests
                
 * [Python][], optional, for data comparison tests

Obtain the source from a repository

    git clone <path-to-repository>/vulture vulture-x.x.x-Source

or by unpacking a source archive

    tar zxvf vulture-x.x.x-Source.tar.gz

Create and enter a build directory

    mkdir vulture-build-linux_$(arch)
    cd vulture-build-linux_$(arch)

Configure the build using cmake, specifying the installation prefix, build type 
and any other compilation options (see below)

    cmake -D CMAKE_INSTALL_PREFIX=/usr/local -D CMAKE_BUILD_TYPE=Release -D WITH_TESTS=ON ../vulture-x.x.x-Source

Build the software

    make

Optionally run the testsuite

    make test

Install the software

    make install

## Developer builds and testing on Linux

The tests can also be run directly using CTest. The software can be built with

    ctest -D ExperimentalBuild
    
and tests run using:

    ctest -R testNames -D ExperimentalTest
    
The `-R` option takes a glob with wildcards of the names of tests to run as its 
argument. Tests can be excluded with `-E testNames`. 

The tests can be run under the valgrind dynamic memory analyser using

    ctest -R testNames -D ExperimentalMemCheck

For example, to run the vulture and gvulture tests with [valgrind][] use

    ctest -R ".*_.*vulture" -D ExperimentalMemCheck

To collect coverage information the software should be configured with build 
type `Pofile`

    cmake -D CMAKE_BUILD_TYPE=Profile ....

and then tests should be run and coverage information collected 

    ctest -R ".*_.*vulture" -D ExperimentalTest
    ctest -R ".*_.*vulture" -D ExperimentalCoverage 

The results will be stored in the `Testing/Temporary` directory.
    
The test results can be submitted to a dashboard using:

    ctest -R ".*_.*vulture" -D ExperimentalSubmit

This is not currently configured.

More verbose output can be obtained using `ctest -V`.

## Building from source on Windows

Requirements: 

 * [Mercurial][], optional, for obtaining source from repository
 
 * [CMake][], required, version >= 2.8.3
 
 * C compiler, , required, with C99 support
 
Obtain the source from a Mercurial repository or by unpacking a source archive 
`vulture-x.x.x-Source.zip` using the windows file manager. 

Create a build directory, for example, `mkdir vulture-build-win32`.

Configure the build using the cmake gui:

 1. Select the source folder.
 2.  Select the build folder created above.
 3.  Run "Configure"
 4.  Set in the GUI
 
     CMAKE_BUILD_TYPE -> Release
     CMAKE_INSTALL_PREFIX -> installation folder, e.g. C:\Program Files\Vulture
     Other options (see below)
 
 5. Run "Configure" again
 6. Run "Generate"

The rest of the compilation process depnds on which compilation suite is been 
used. For the Mingw compiler:

 1. Open a DOS terminal in the build folder.
 2. Compile the software

    mingw32-make

 3. Run the testsuite

    mingw32-make test

 4. Install the software

    ming32-make install


## Compilation options

The following complation options are supported.

### USE_SCALED_FIELD=ON/OFF (default: OFF)

Use scaled electric and magnetic fields. Each field component is scaled by its 
corresponding grid edge length. This process is transparent to the user. It 
allows for slightly faster code by redcung the operation count in the main field 
update loops. This option is not compatiable with `USE_INDEXED_MEDIA`.

### USE_AVERAGED_MEDIA=ON/OFF (default: OFF)

Average simple media parameters across material interfaces. This is more 
accutate than the direct per voxel setting of the material parameters. This 
options is not compatiable with `USE_INDEXED_MEDIA`.

### USE_INDEXED_MEDIA=ON/OFF (default: OFF)

Use integer indices to hold the medium type a each field point rather than real 
arrays of update coefficients. This can significantly reduce the memory 
requirements at the expense of a slight decrease in performance. The size the 
interger index is dset usiing the type MediumIndex and macro `MAX_MEDIA UINT_MAX` 
in header file medium.h. By default this is an unsigned int but could be an 
unsigned char to minimise the memory requirement. This option is not compatible 
with `USE_SCALED_FIELD` or `USE_AVERAGED_MEDIA`.

### WITH_OPENMP=ON/OFF (default: OFF)

Use [OpenMP][] directives to parallelise the main field update loops. The number of 
threads to use is set by the environment variable `OMP_NUM_THREADS`. On Linux this 
can be defined from the shell using

    export OMP_NUM_THREADS=2
    
for bash/sh or

    setenv OMP_NUM_THREADS 2
    
for csh/tcsh. Under Windows it can be set from the start Menu.

### WITH_SIBC=ON/OFF (default: OFF)

Build with surface impedance boundary condition support.

### WITH_TESTS=ON/OFF (default: OFF)

Build with support for running the test suite.
  
### GNUPLOT_TESTS=ON/OFF (default: OFF)

Include testsuite components that require gnuplot if it is available.

### PYTHON_TESTS=On/OFF (default: OFF)

Include testsuite components that require python if it is available.

###OCTAVE_TESTS=ON/OFF (default: OFF)

Include testsuite components that require GNU Octave if it is available. 

### PROCESSING_TESTS=ON/OFF (default: OFF)

Include testsuite components that require the AEG post-processing package if it 
is available. 

### BIG_TESTS=ON/OFF (default: OFF)

Include tests that take a long time to run.

### WITH_LATEX=ON/OFF (default: OFF)

Build the LaTeX user manual. 



[git]: https://git-scm.com
[CMake]: http://www.cmake.org
[Octave]: http://www.octave.org
[MATLAB]: http://www.mathworks.com/MATLAB
[Python]: http://www.python.org
[valgrind]: http://valgrind.org
[OpenMP]: http://openmp.org/wp
