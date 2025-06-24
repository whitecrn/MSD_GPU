#ifndef MSD_H_
#define MSD_H_

#include "error.cuh"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <cuda_runtime.h>
#include "type.h"

using namespace std;


void read_dump(vector<atom_old>& R_old,int &atom_number,int &step,double &a,double &b,double &c);
void __global__ moving_xyz(atom *,const double,const double,const double,const int,const int,const int);
void __global__ MSD_cal(atom *,double *,const int,const int,const int);
void MSD_write(double *);
int show_information();

#endif