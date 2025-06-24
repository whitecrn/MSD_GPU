#ifndef TYPE_H_
#define TYPE_H_

extern double potim;
extern int atom_id;
extern int step;

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

#endif