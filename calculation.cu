#include "data.h"

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                         MSD_CALCULATION                              //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
void __global__ MSD_cal
(atom *d_R,double *d_msd,const int step,const int atom_number,const int number)
{   
    const int n1=blockIdx.x*blockDim.x+threadIdx.x;
    const int n2=blockIdx.y*blockDim.y+threadIdx.y;
    const int delta_t=static_cast<int>(round(0.5*step));
    const int t=delta_t-1;
    if ((n1<t) && (n2<delta_t))
    {   
        for (int i=0;i<number;i++)
        {
            double dx=d_R[n1*atom_number+i].x-d_R[(n1+n2)*atom_number+i].x;
            dx*=dx;
            double dy=d_R[n1*atom_number+i].y-d_R[(n1+n2)*atom_number+i].y;
            dy*=dy;
            double dz=d_R[n1*atom_number+i].z-d_R[(n1+n2)*atom_number+i].z;
            dz*=dz;
            double dr=dx+dy+dz;
            atomicAdd(&d_msd[n2],dr);
        }
    }
}