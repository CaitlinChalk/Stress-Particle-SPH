# Stress-Particle-SPH
Code and documentation for the Stress-Particle SPH method, developed in the PhD project of Caitlin Chalk

Contained are the raw Fortran files (.f90), with example problems and documentation.
NOTE: The .f90 files are included separately for each problem. 
      This is because the code is still under development, and needs to be tested for consistency amongst the problems.
The documentation includes the PhD thesis for which the work was developed, and descriptions of the input files.
Further documentation is to be added.

File information:
Fortran files: .f90
Input files: input.txt, filename.dat, filename.pts
Output files: output.txt, filename.chk, .csv files (for ParaView visualisation), filename.post.res, filename.post.msh (for GiD visualisation)
Make file: makefile (for compilation purposes)

Instructions to run on a linux machine:

1) Ensure the .f90 and makefile files are in the same directory
2) Clean the directory by typing 'make clean' in the terminal.
3) Compile the code by making 'make' in the terminal. An executable named 'sph' will appear
4) Copy the execute sph file to the same directory as the three input files
5) Run the code by typing './sph' in the terminal.
