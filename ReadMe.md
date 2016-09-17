![](https://bitbucket.org/uoyaeg/vulture/wiki/aperture.jpg )

[TOC]

# Vulture: An Open Source FDTD Solver For Electromagnetic Simulations

The Applied Electromagnetics Group ([AEG][]) [FDTD][] solver, Vulture, is an 
[Open Source][] non-uniform structured mesh [FDTD][] code for electromagnetic 
simulations. It was developed in the in the [Department of Electronics][] at the 
[University of York][] for research in electromagnetic compatibility ([EMC][]), 
computational electromagnetics ([CEM][]) and [bioelectromagnetics][].

## Code Features

The code currently has the following features:

* Non-uniform mesh allowing for uniform cubic and uniform cuboid special cases.

* External mesh surfaces that can independently be perfect electric conductor 
  (PEC), perfect magnetic conductor (PMC), perfectly match layer ([PML][]), analytic 
  Mur absorbing boundary condition (ABC) or periodic boundary conditions.
 
* A uniaxial perfectly matched layer (UPML) implementation that can terminate 
  arbitrary inhomogeneous media.

* Gaussian pulse, compact pulse, ramped sinusoid, differentiated pulse and user 
  defined waveforms.

* Distributed hard and soft electric and magnetic field, current density, 
  current and ideal voltage sources.

* Lumped resistive voltage and current sources.

* Internal PEC surfaces.

* Simple isotropic media with frequency independent permittivity, conductivity 
  and (real) permeability.

* Arbitrary electrically dispersive media using a generalised multi-pole Debye 
  dispersion relationship.

* A total-field-scattered-field (TFSF) plane-wave source, also known as a 
  Huygen's surface source, for multiple plane-wave excitation. The implementation 
  supports partial Huygen's surfaces and has grid dispersion optimisations for 
  uniform cubic meshes.

* Binary and ASCII format field observers.
 
The mesh format is a simple ASCII file that can be written manually or for more 
complex models [AEG Mesher][] can export *Vulture* format meshes.

## Requirements

The code is written in standard [C99][]. Additional requirements are:

1. (Mandatory) To compile and install the code the [CMake][] software build tool 
   is needed.

2. (Optional) To help with development or as an alternative way to download the 
   source a client for the [Mercurial][] Version Control System is required.

To run the test-suite the following are also required:

1. (Mandatory) GNU [Octave][] or [MATLAB][].

2. (Optional) [gnuplot][].

3. (Optional) [Python][].

4. (Optional) The [AEG][] time-domain post-processing tools. 

The code has been primarily developed on Linux platforms, but it should build 
and run on both Linux and Windows systems.

## Documentation

Installation instructions are contained in the file [Install.txt][] in the 
source distribution. There is also a [LaTeX][] user manual and tutorial in the 
doc directory of the source distribution that builds into [UserManual.pdf][].

## Bugs and support

The code is still under development and no doubt will contain many bugs. Known 
significant bugs are listed in the file doc/[Bugs.txt][] in the source code. 

Please report bugs using the bitbucket issue tracker at
<https://bitbucket.org/uoyaeg/vulture/issues> or by email to <ian.flintoft@york.ac.uk>.

For general guidance on how to write a good bug report see, for example:

* <http://www.chiark.greenend.org.uk/~sgtatham/bugs.html>
* <http://noverse.com/blog/2012/06/how-to-write-a-good-bug-report>
* <http://www.softwaretestinghelp.com/how-to-write-good-bug-report>

Some of the tips in <http://www.catb.org/esr/faqs/smart-questions.html> are also 
relevant to reporting bugs.

There is a Wiki on the bitbucket [project page](https://bitbucket.org/uoyaeg/vulture/wiki).

## How to contribute

We welcome any contributions to the development of the code, including:

* Fixing bugs.

* Interesting examples that can be used for test-cases.

* Improving the user documentation.

* Items in the to-do list in the file [ToDo.txt][].

Please contact [Dr Ian Flintoft], <ian.flintoft@york.ac.uk>, if you are 
interested in helping with these or any other aspect of development.

## Licence

The code is licensed under the [GNU Public Licence, version 3](http://www.gnu.org/copyleft/gpl.html). 

## Developers

[Dr Ian Flintoft](http://www.elec.york.ac.uk/staff/idf1.html) : <ian.flintoft@york.ac.uk>

Mr Samuel Bourke : <sab544@york.ac.uk>

## Contacts

[Dr Ian Flintoft](http://www.elec.york.ac.uk/staff/idf1.html) : <ian.flintoft@york.ac.uk>

[Dr John Dawson](http://www.elec.york.ac.uk/staff/jfd1.html) : <john.dawson@york.ac.uk>

## Credits

Vulture was inspired by a number of earlier electromagnetic simulation codes 
developed in the Applied Electromagnetics Group in the [Department of 
Electronics][] at the [University of York][]. The input mesh format for Vulture 
is an evolution of the format used by the “Hawk” Transmission Line Matrix 
([TLM][]) and “Falcon” [FDTD][] solvers by [Dr Stuart Porter][] and [Dr John 
Dawson][]. The binary output format was also originally developed by [Dr Stuart 
Porter][] and [Dr John Dawson][] for these two codes.

## Publications using Vulture

[Xia2012]: http://dx.doi.org/10.1109/EMCEurope.2012.6396718

([Xia2012]) R. Xia, J. F. Dawson, I. D. Flintoft, A. C. Marvin and S. J. Porter, “Use of a 
genetic algorithm in modelling small structures in airframe”, EMC Europe 2012, 11th 
International Symposium on EMC, Rome, Italy, 17-21 September, 2012.

[Xia2011]: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6078646&isnumber=6078493

([Xia2011]) R. Xia, J. F. Dawson, I. D. Flintoft, A. C. Marvin, S. J. Porter and I. Marschke, 
“Building macro-models for small structures on aircraft”, EMC Europe 2011, 10th International 
Symposium on EMC, York, UK, 26-30 September, 2011. pp. 575-580. 

## Related links

* The free [MEEP](http://ab-initio.mit.edu/wiki/index.php/Meep) devloped at MIT.

* [Understanding the FDTD Method](http://www.eecs.wsu.edu/~schneidj/ufdtd) by John B. Schneider. 


[PML]:                       http://en.wikipedia.org/wiki/Perfectly_matched_layer
[Open Source]:               http://opensource.org
[LaTeX]:                     http://www.latex-project.org
[TLM]:                       http://en.wikipedia.org/wiki/Transmission-line_matrix_method
[FDTD]:                      http://en.wikipedia.org/wiki/Finite-difference_time-domain_method
[gnuplot]:                   http://www.gnuplot.info
[Python]:                    https://www.python.org
[Octave]:                    http://www.gnu.org/software/octave
[MATLAB]:                    http://www.mathworks.co.uk/products/matlab
[C99]:                       http://en.wikipedia.org/wiki/C99
[Mercurial]:                 http://mercurial.selenic.com
[CMake]:                     http://www.cmake.org
[AEG Mesher]:                https://bitbucket.org/uoyaeg/aegmesher

[University of York]:        http://www.york.ac.uk
[Department of Electronics]: http://www.elec.york.ac.uk
[AEG]:                       http://www.elec.york.ac.uk/research/physLayer/appliedEM.html
[Dr Ian Flintoft]:           http://www.elec.york.ac.uk/staff/idf1.html
[Dr John Dawson]:            http://www.elec.york.ac.uk/staff/jfd1.html
[Dr Stuart Porter]:          http://www.elec.york.ac.uk/staff/sjp1.html
[EMC]:                       http://www.elec.york.ac.uk/research/physLayer/appliedEM/emc.html
[CEM]:                       http://www.elec.york.ac.uk/research/physLayer/appliedEM/numerical.html
[bioelectromagnetics]:       http://www.elec.york.ac.uk/research/physLayer/appliedEM/bio.html

[Install.md]:                https://bitbucket.org/uoyaeg/vulture/raw/tip/Install.md
[UserManual.pdf]:            https://bitbucket.org/uoyaeg/vulture/wiki/UserManual.pdf
[ToDo.md]:                   https://bitbucket.org/uoyaeg/vulture/raw/tip/doc/ToDo.md
[Bugs.md]:                   https://bitbucket.org/uoyaeg/vulture/raw/tip/doc/Bugs.md