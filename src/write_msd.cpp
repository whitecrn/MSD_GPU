#include "../include/write_msd.h"
#include "../include/type.h"

void MSD_write(double *h_msd_ave)
{
    ofstream write_file;
    write_file.open("data/msd.dat");
    write_file.precision(9);
    write_file.setf(ios::fixed, ios::floatfield);
    if (!write_file.is_open())
    {
        cout << "Can't open msd.dat file." << endl;
        exit(1);
    }

    for (int t=0;t<round(0.5*step);t++)
    {
        //cout.precision(9);
        //cout.setf(ios::fixed,ios::floatfield);
        write_file << double(potim*double(t)) << "  " << h_msd_ave[t]  << endl;
    }
    write_file.close();
}