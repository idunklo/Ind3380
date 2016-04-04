To complete oblig 1 you have to first write a serial program
by filling out the serial_main.c file in the serial/
directory and finally write a parallel program by filling
out the parallel_main.c file in the parallel/ directory.

The two programs shall perform an iso diffusion denoising on
a JPEG image. The two programs shall also accept the
following parameters (in order)

$ ./program number_of_iterations kappa_value infile outfile

When you are in the serial/ or parallel/ directories you may
compile the programs by typing "make" in the terminal.

When you are ready to deliver you enter the base directory
and type "make delivery" in the terminal. A tarball with
your updated files will then be created ready for delivery.

##########################################################
Added to readme by me:

To run program for parallell:
$ mpirun -np 4 ./program number of iterations kappa_value infile outfile

Have cooperated with Mathias Kirkerød. The program doesn't work with how you can't work with the same process for the 0-process as for the others, but I can't figure out how to fix it.
