

Vulture: An Open Source FDTD Solver For Electromagnetic Simulations

The Applied Electromagnetics Group (AEG) FDTD solver, Vulture, is an Open_Source
non-uniform structured mesh FDTD code for electromagnetic simulations. It was
developed in the in the Department_of_Electronics at the University_of_York for
research in electromagnetic_compatibility_(EMC), computational_electromagnetics_
(CEM) and computational_dosimetry_and_bioelectromagnetics.

Code Features

The code currently has the following features:

* Non-uniform mesh allowing for uniform cubic and uniform cuboid special cases.
* External mesh surfaces that can independently be perfect electric conductor
  (PEC), perfect magnetic conductor (PMC), perfectly match layer (PML), analytic
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
  Huygen's surface source, for multiple plane-wave excitation. The
  implementation supports partial Huygen's surfaces and has grid dispersion
  optimisations for uniform cubic meshes.
* Binary and ASCII format field observers.

The mesh format is a simple ASCII file that can be written manually or for more
complex models AEG_Mesher can export Vulture format meshes.

Requirements

The code is written in standard C99. Additional requirements are:

  1. (Mandatory) To compile and install the code the CMake software build tool
     is needed.
  2. (Optional) To help with development or as an alternative way to download
     the source a client for the Mercurial Version Control System is required.

To run the test-suite the following are also required:

  1. (Mandatory) GNU Octave or MATLAB.
  2. (Optional) gnuplot.
  3. (Optional) Python.
  4. (Optional) The AEG time-domain post-processing tools.

The code has been primarily developed on Linux platforms, but it should build
and run on both Linux and Windows systems.

Documentation

Installation instructions are contained in the file Install.txt in the source
distribution. There is also a LaTeX user manual and tutorial in the doc
directory of the source distribution that builds into a PDF_file.

Bugs and support

The code is still under development and no doubt will contain many bugs. Known
significant bugs are listed in the file doc/Bugs.txt in the source code.
Please report bugs using the bitbucket issue tracker at https://bitbucket.org/
uoyaeg/vulture/issues or by email to ian.flintoft@york.ac.uk.
For general guidance on how to write a good bug report see, for example:

* http://www.chiark.greenend.org.uk/~sgtatham/bugs.html
* http://noverse.com/blog/2012/06/how-to-write-a-good-bug-report
* http://www.softwaretestinghelp.com/how-to-write-good-bug-report

Some of the tips in http://www.catb.org/esr/faqs/smart-questions.html are also
relevant to reporting bugs.
There is a Wiki on the bitbucket project_page.

How to contribute

We welcome any contributions to the development of the code, including:

* Fixing bugs.
* Interesting examples that can be used for test-cases.
* Improving the user documentation.
* Items in the to-do list in the file ToDo.txt.

Please contact [Dr Ian Flintoft], ian.flintoft@york.ac.uk, if you are interested
in helping with these or any other aspect of development.

Licence

The code is licensed under the GNU_Public_Licence,_version_3.

Developers

Dr_Ian_Flintoft : ian.flintoft@york.ac.uk
Mr Samuel Bourke : sab544@york.ac.uk

Contacts

Dr_Ian_Flintoft : ian.flintoft@york.ac.uk
Dr_John_Dawson : john.dawson@york.ac.uk

Credits

Vulture was inspired by a number of earlier electromagnetic simulation codes
developed in the Applied Electromagnetics Group in the Department_of_Electronics
at the University_of_York. The input mesh format for Vulture is an evolution of
the format used by the “Hawk” Transmission Line Matrix (TLM) and
“Falcon” FDTD solvers by Dr_Stuart_Porter and Dr_John_Dawson. The binary
output format was also originally developed by Dr_Stuart_Porter and Dr_John
Dawson for these two codes.

Related links


* The free MEEP devloped at MIT.
* Understanding_the_FDTD_Method by John B. Schneider.

