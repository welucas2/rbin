# README #

rbin is a Fortran 90 module with which you can quickly read and store data for use from an sphNG binary file.

### How to use ###

The code only needs a modern Fortran compiler. Any recent version of gfortran or ifort should do.

The file should be opened outside the module in your own code with logical unit 10, as an 'unformatted' F77 file.

The logical variable 'contiguous' should be set in your own code before reading when working with files using MPI blocks. When .TRUE., the module will read the MPI blocks sequentially, storing the arrays from each block in single contiguous arrays. As such, the entire dataset can be accessed without further read operations. If the particle number is large enough that your machine is not able to read the entire file, set contiguous = .FALSE. In this case, control will pass back to you after each block is read. In this case you will be reading one block, performing your calculations, then reading the next block, and so on.

You will see at the top of the module several parameters 'igradh', 'imhd', 'iexf' and 'imigrate'. These should be set to 0 or 1 as needed for your data. Note though that the code may have issues with these as they aren't particularly well implemented yet!

### For the future ###

Re-use arrays - rather than reallocating for each new block, set them equal in size to the largest block and vary maximum particle index instead. This may help with memory issues on OSX. Opening multiple files in succession could still be problematic.

Move away from a static set of data towards something more more object oriented. This would be good for e.g. comparing multiple files side-by-side.

### Who do I talk to? ###

Email wel2@st-andrews.ac.uk for help.