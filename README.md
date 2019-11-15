# vibAnalysis.exe

The program `vibAnalysis.exe` carries out vibration analysis given as set of atomic masses, atomic Cartesian coordinates, and a Hessian matrix. This program calculates vibrational frequencies (reported in cm^-1^) and normalize un-mass-weighted normal displacement vectors. By default, the solution of the vibrational analysis problem includes projection of overal translational and rotational degrees of freedom. Optionally, the user may indicate 1 or more frozen atomic centers. In such a case, the program applies the correct projections on the Hessian matrix to remove atomic center translations from the subsequent normal mode formation.

## Usage

To use this program, one must have produced a Gaussian matrix file. For help doing that, one should consult the Gaussian User Guide [available online](http://gaussian.com/interfacing/). An example command for a Gaussian matrix file named _test.mat_ is

````shell
./vibAnalysis.exe test.mat
````

To carry out vibration analysis with frozen atoms, one gives a list of atomic centers frozen (by their center number) after the matrix file name. For example, to carry out a vibrational analysis on the system defined in matrix file _test.mat_ with atomic centers 2, 5, and 12 frozen one would use the command

````shell
./vibAnalysis.exe test.mat 2 5 12
````

## Compiling the Program

This program depends on the MQCPack library. This library can be found on the [MQC site](https://github.com/MQCPack). The makefile assumes the user has the standard $mqcroot variable set in their environment. As an alternative set-up, the makefile can be modified accordingly.

The makefile also assumes that the `pgfortran` command is available. Using other fortran 2003 compilers should be possible, but appropriate modifications to the makefile will be required.

With the MQCPack library built and the PGI compiler available, vibAnalysis.exe is build by running `make`.

## Sub-Directories

There is a sub-directory with test jobs and outputs named GTests. Gaussian input files and formatted checkpoint files are both available in that directory. One can run the input files to generate Gaussian outputs and matrix files, or one can use the Gaussian utility `formchk` to build matrix files from the formatted checkpoint files. Information on that process is available on the [Gaussian website](http://gaussian.com/formchk/).

# Author

This program was written by Hrant P. Hratchian, University of California, Merced. Questions and bug reports can be sent by email to hhratchian@ucmerced.edu.

# Version History

The latest version of this program was uploaded to GitHub in November, 2019.