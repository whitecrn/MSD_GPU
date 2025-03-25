#include "error.cuh"
#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <cuda_runtime.h>
using namespace std;

extern double potim;
extern int atom_id;
extern int step;
extern string file_name;

struct atom_old {
    int step;  
    int indices;  
    int atom_type; 
    double x, y, z; 

    atom_old(int s, int idx, int t, double x_, double y_, double z_)
        : step(s), indices(idx), atom_type(t), x(x_), y(y_), z(z_) {}

    bool operator<(const atom_old& other) const 
    {
        if (step != other.step) 
        {
            return step < other.step;
        }
        return indices < other.indices; 
    }
};

struct atom
{
    int atom_type;
    double x,y,z;
};

const int TILE_DIM=32;

void read_dump(vector<atom_old>& R_old,int &atom_number,int &step,double &a,double &b,double &c);
void __global__ moving_xyz(atom *,const double,const double,const double,const int,const int,const int);
void __global__ MSD_cal(atom *,double *,const int,const int,const int);
void MSD_write(double *);