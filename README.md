# WOSS
WOSS (wavelet optimization by stochastic switch) is an experimental algorithm to build or improve wavelets in reservoir engineering. The algorithm assumes an initial wavelet (ricker is ok) well log individual files and a seismic cube (currently reading numpy binary [npy]).

## How does it work
It reads a wavelet file (it can be a basic Ricker wavelet), a seismic cube (currently in 3D npy format), and a list of well log files. For every iteration it will randomly switch a value in the wavelet and convolves the wells. If it gets better the switch is accepted, otherwise it won't be considered. In the end it will output a plot with the comparison between the original and optimized wavelet as well as the initial and final correlation (quasi-correlation is being used). I didn't implement a file writer for the final wavelet but it should be pretty easy to do so (just save wathever it is on final_wavelet variable).

## What's a wavelet file
This is an example of a wavelet file:

-4	-150410
-3	-188835
-2	-109970
-1	 121088
 0	 401148
 1	 545489
 2	 451922
 3	 191310
 4	-54597.1
 
 First column is just the step (ignore it but put it there), the second the actual wavelet signal. The wavelet size is actually the number of points until you get to 0 in the wavelet step. In our example is 4.
 
 ## What is a well log file
 A well log file is a file with simple X,Y,Z,AI columned information such as this:
 
      326.66667       144.66667         0.25000      9278.43000
      326.66667       144.66667         1.25000      9278.24000
      326.66667       144.66667         2.25000      9278.05000
      326.66667       144.66667         3.25000      9277.86000
      326.66667       144.66667         4.25000      9277.66000
      326.33333       144.66667         5.25000      9277.47000
      326.33333       144.66667         6.25000      9277.28000
      326.33333       144.66667         7.25000      9277.09000
      326.33333       144.66667         8.25000      9276.89000
      326.33333       144.66667         9.25000     13694.50000
      326.33333       144.66667        10.25000     14945.32000
      
Please notice that the well should come in regular GRID coordinates so x=326 is actually cell row=326 in the seismic. Not cell size or first coordinate is being considered for this version of the algorithm. All of this is important because you'll probably have to transform your data to use this algorithm.

## What is a seismic file
The seismic file is binary and I could build an analogue in Python using the following code:

import numpy as np
g = np.zeros((100,100,30),dtype='float32') # So building a 3D matrix of type float32 with number of nodes "x,y,x" of (100,100,30).
np.save('zero_matrix.npy,g)                # I'm just saving a 3D cube with zeros inside it to a binary file. You should actually populate it with your seismic values.

