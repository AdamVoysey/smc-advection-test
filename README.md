
	Guide notes to run this SMC grid propagation model
                    Jian-Guo Li   18 July 2019

The propagation test model is has been set up in order to trial the SMC
wave propagation module on GPUs or other HPC architectures.

The source codes are stored in CodSrc and they consist of 4 modules:
WRefractnHyb.f90:	The main program and some auxilary subroutines.
W3PSMCMD.f90:		Main SMC grid module in WW3 but simplified for the test.
W3DispMD.f90:		Original WW3 module for calculation of wave parameters.
Constants.f90:		Data structure and particular grid and model settings.

The program has presently been compiled and run on the Met Office xcs machine,
using the Makefile in CodSrc.  You may need to modify the Makefile if you are on a
different machine.

The SMC grid used in this test is the 3-6-12-25 km four level SMC grid for the
Mediterranean Sea and the cell and face array files are stored in DatGMC/.

The jobscript smcrHybMed30N is for job queue on xcs. It uses 30 nodes in a
hybrid configuration of 4 ranks per node and 9 threads per rank, resulting in a
total of 120 ranks and 1080 PEs.  You may try a 12 ranks per node and 3 threads
per rank option with this job script by simplying swapping the corresponding
lines in the script, namely the OMP_NUM_THREADS and aprun lines.

The jobscript smcrHybMed01N is for job queue on one single xcs node.  There are
3 tested options for the 36 cores with 3, 9 and 18 threads or 12, 4 and 2 ranks.
You may select them by commenting off other option lines as in the 30 nodes one.

The model results are stored in OutDat as text files.  They are total wave
height at all sea points every 2 hr model time.  These data files could be
visualised with Python programs stored in PyPlot.

Execute python smed36125grids.py on xcs will generate the Med36125 grid plot as
a ps file.  Most importantly, it will produce a mapping or projection file and
store it in DatGMC.  This is a binary file for python to read and should be
created on your local machine just in case the binary format is different.

Then try python smed36125props.py to draw the SWH plots from the test results.
If you just run the model on your machine to generate the results, please use
the linux command line 
PropHybr/OutDat/  ls C* > cfiles.txt
to produce a list of files for the python program to draw. You may select some 
of output files for plots rather than all of them.

Before you run the model and the python program, please go through them to
change the paths to your own ones.  The CodSrc/Constants.f90 has some paths for
the cell and face array files which you need to modify.  The job script also
has some input and output directories, which you need to point to yours.  And
finally the two python programs also contain paths to input and output files.

One more thing you need to note is that there are a few python library files
stored in PyPlot.  You may eigher move them into your python library folder or
simply run the python program in the directory of PyPlot.

Any queries or suggestions to improve the model please feel free to email me at 
Jian-Guo.Li@metoffice.gov.uk

Good luck and happy testing the SMC grid propagation model. 

