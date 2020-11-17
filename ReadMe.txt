Author: 
Shashank Dattathri
UG 4th Year
#13799

My project uses the following packages: numpy, scipy, matplotlib, and astropy. The first three are most probably installed on most systems used for scientific computing. 
To install astropy, one simply needs to type in "pip install astropy" in the command line. 

The folder Data contains the data files used in my project. They are all of the .fit format. Note that these are not the original published data; I significantly cropped
them and changed their format so that my submission folder is of a reasonable size. I did not further reduce them to regular data tables because reading and writing FITS
files was an essential part of my project and it was non-trivial to learn. 

The folder Codes contains 8 python files. The file azimuthal_average.py is used in each of the other calculations, and does not run by itself. The other six files that 
contain my calculations are Planck.py, WISE.py, TB_vs_vlsr.py, Col_density_plot.py, Multipase_ISM.py, and Solenoidal_fraction.py. The file main.py executes each of the
six files one by one. This is the easiest method to test my code. Alternatively, one can run each of the six files individually. 

The full description of what each of the files calculates is given in my report. 