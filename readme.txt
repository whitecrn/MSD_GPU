a script calculate average mean square displacement.

how to use:
1. nvcc command is available. 

2. make

3. ./msd_gpu -i <atom type for calculate> -f <input file> -p <potim> -h for help -t <type of coordinate>

atom id: atom id used to calculated
input file: input file name
potim: unit ps
type of coordinate: Cartesian or Fraction. 0 for Fraction, other number for Cartesian.


default: 
atom id: 1
input file: dump.atom
potim: 0.1 ps
type of coordinate: Cartesian
