a script calculating average mean square displacement.

how to use:
1. ensure that you have installed CUDA and nvcc command is available. 

2. make

3. ./msd_gpu -i <atom id> -f <input file> -p <potim> -h for help

atom id: atom id used to calculated
input file: input file name
ppotim: unit ps

default: 
atom id: 1
input file: dump.atom
potim: 0.1 ps
