#include "../include/boundary.h"
#include "../include/type.h"

void __global__ moving_xyz
(atom *d_R,const double a,const double b,const double c,const int atom_number,const int step,const int atom_id)
{
    const int n=blockIdx.x*blockDim.x+threadIdx.x;
    if ((n<atom_number) && (d_R[n].atom_type==atom_id))
    {
        for (int i=0;i<step-1;i++)
        {
            for (int j=0;j<100;j++)
            {
                double dx=d_R[i*atom_number+n].x-d_R[(i+1)*atom_number+n].x;
                double dy=d_R[i*atom_number+n].y-d_R[(i+1)*atom_number+n].y;
                double dz=d_R[i*atom_number+n].z-d_R[(i+1)*atom_number+n].z; 
                if (dx>=(0.5*a)) d_R[(i+1)*atom_number+n].x+=a;
                else if (dx<=(-0.5*a)) d_R[(i+1)*atom_number+n].x-=a;

                if (dy>=(0.5*b)) d_R[(i+1)*atom_number+n].y+=b;
                else if (dy<=(-0.5*b)) d_R[(i+1)*atom_number+n].y-=b;

                if (dz>=(0.5*c)) d_R[(i+1)*atom_number+n].z+=c;
                else if (dz<=(-0.5*c)) d_R[(i+1)*atom_number+n].z-=c;
            }
        }
    } 
}