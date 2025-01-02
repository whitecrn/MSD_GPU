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
using namespace std;

struct atom_old {
    int step;  
    int indices;  
    int atom_type; 
    double x, y, z; 

    atom_old(int s, int idx, int t, double x_, double y_, double z_)
        : step(s), indices(idx), atom_type(t), x(x_), y(y_), z(z_) {}

    bool operator<(const atom_old& other) const {
        if (step != other.step) {
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


double potim=0.1;
int atom_id=0;
const int TILE_DIM=32;
int step=0;
string file_name="dump.atom";
bool h_flag=false;

void read_dump(vector<atom_old>& R_old,int &atom_number,int &step,double &a,double &b,double &c);
void __global__ moving_xyz(atom *,const double,const double,const double,const int,const int,const int);
void __global__ MSD_cal(atom *,double *,const int,const int,const int);
void MSD_write(double *);

int main(int argc,char *argv[])
{   
    int device_id=-1;
    cudaError_t err=cudaGetDevice(&device_id);
    if (err != cudaSuccess) 
    {
        cout << "CUDA error: Unable to get current device. " << cudaGetErrorString(err) << endl;
        return 1;
    }
    cudaDeviceProp device_prop;
    err=cudaGetDeviceProperties(&device_prop,device_id);
    if (err != cudaSuccess) 
    {
        cout << "CUDA error: Unable to get device properties. " << cudaGetErrorString(err) << endl;
        return;
    }

    cout << "Device " << device_id << ": " << device_prop.name << endl;
    cout << "Compute capability: " << device_prop.major << "." << device_prop.minor << endl;
    int opt;
    while ((opt=getopt(argc,argv,"i:p:f:h"))!=-1)
    {
        switch (opt)
        {
            case 'i':
                atom_id=atoi(optarg)-1;
                break;
            case 'p':
                potim=atof(optarg);
                break;
            case 'f':
                file_name=optarg;
                break;
            case 'h':
                h_flag=true;
                break;
            default:
                cout << "Usage: " << argv[0] << " " << "-i <Atom ID> -p <Potim> -f <File Name> -h <Help>" << endl;
                return 1;
        }
    }

    if (h_flag)
    {
        cout << "Usage: " << argv[0] << "-i <Atom ID> -p <Potim> -f <File Name> -h <Help>" << endl;
        cout << "-i <Atom ID>: Set Atom ID(int) which will use to calculate MSD" << endl;
        cout << "-p <Potim>: Set Potim(double) of your system, unit: ps" << endl;
        cout << "-f <File Name>: Set name of your input dump file(string)" << endl;
    }
    cout << "Potim: " << potim << endl;
    cout << "Atom_ID: " << atom_id+1 << endl;
    cout << "File Name: " << file_name << endl;

    double a,b,c;
    long N,N_atom;
    int atom_number=0;
    vector<atom_old> R_old;

    clock_t start1=clock();
    read_dump(R_old,atom_number,step,a,b,c); // read dump file.
    int Max=0;
    for (int i=0;i<atom_number;i++)
    {
        if (R_old[i].atom_type>Max)
            Max=R_old[i].atom_type;
    }
    clock_t end1=clock();

    double duration=double(end1-start1)/CLOCKS_PER_SEC;
    cout << "Finished reading dump file, time: " << duration << " s" << endl;


    N=R_old.size(); //R and R_old size.
    N_atom=N*sizeof(atom);
    atom *R=new atom[N];
    int number[Max+1];
    for (int i=0;i<Max;i++) number[i]=0;
    for (int i=0;i<N;i++) //copy R_old to R.
    {
        R[i].atom_type=R_old[i].atom_type;
        if (i<atom_number) 
        {
            number[R[i].atom_type]+=1;
        }
        R[i].x=R_old[i].x;
        R[i].y=R_old[i].y;
        R[i].z=R_old[i].z;
    }
    int deltaT=static_cast<int>(round(0.5*step)); //delta_t size.
    int msd_size=sizeof(double)*deltaT;
    double *h_msd=new double[deltaT];
    double *h_msd_ave=new double[deltaT];
    double *d_msd;
    atom *d_R;

    memset(h_msd,0.0,deltaT*sizeof(double)); //initial h_msd.
    memset(h_msd_ave,0.0,deltaT*sizeof(double)); //initial h_msd_ave.

    CHECK(cudaMalloc(&d_msd,msd_size)); //allocate d_msd with round(0.5*step)*sizeof(double).
    CHECK(cudaMalloc(&d_R,N_atom)); //allocate d_R with R size*sizeof(atom).

    CHECK(cudaMemcpy(d_R,R,N_atom,cudaMemcpyHostToDevice)); // from R to d_R
    CHECK(cudaMemcpy(d_msd,h_msd,msd_size,cudaMemcpyHostToDevice)); // from h_msd to d_msd.
//moving xyz:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    cudaEvent_t start,stop;
    CHECK(cudaEventCreate(&start));
    CHECK(cudaEventCreate(&stop));
    CHECK(cudaEventRecord(start));
    cudaEventQuery(start);

    moving_xyz<<<(N+128-1)/128,128>>>(d_R,a,b,c,atom_number,step,atom_id); //moving xyz.
    cudaDeviceSynchronize();

    CHECK(cudaEventRecord(stop));
    CHECK(cudaEventSynchronize(stop));
    float elapsed_time;
    CHECK(cudaEventElapsedTime(&elapsed_time,start,stop));
    cout << "End of moving xyz, time: " << elapsed_time/1000 << " s." << endl;
    CHECK(cudaEventDestroy(start));
    CHECK(cudaEventDestroy(stop));

//calculating MSD:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHECK(cudaEventCreate(&start));
    CHECK(cudaEventCreate(&stop));
    CHECK(cudaEventRecord(start));
    cudaEventQuery(start);

    const int grid_size_x=(deltaT+TILE_DIM-1)/(TILE_DIM);
    const int grid_size_y=grid_size_x;
    const dim3 block_size(TILE_DIM,TILE_DIM);
    const dim3 grid_size(grid_size_x,grid_size_y);

    MSD_cal<<<grid_size,block_size>>>(d_R,d_msd,step,atom_number,number[atom_id]); //MSD calculation.
    cudaDeviceSynchronize();

    CHECK(cudaEventRecord(stop));
    CHECK(cudaEventSynchronize(stop));
    CHECK(cudaEventElapsedTime(&elapsed_time,start,stop));
    cout << "End of calculating, time: " << elapsed_time/1000 << " s." << endl;
    CHECK(cudaEventDestroy(start));
    CHECK(cudaEventDestroy(stop));
//for average:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHECK(cudaMemcpy(h_msd,d_msd,msd_size,cudaMemcpyDeviceToHost)); // move d_msd to h_msd.
    int n=number[atom_id]*(static_cast<int>(0.5*step)-1);
    for (int i=0;i<deltaT;i++) 
    {
        h_msd_ave[i]=h_msd[i];
        h_msd_ave[i]/=n;
    }
//writing MSD:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    MSD_write(h_msd_ave);
// free memory:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHECK(cudaFree(d_msd));
    CHECK(cudaFree(d_R));
    delete [] h_msd;
    delete [] h_msd_ave;
    delete [] R;
    return 0;
}
//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                          dump read                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
void read_dump
(vector<atom_old>& R_old,int &atom_number,int &step,double &a,double &b,double &c)
{    
    double x_b[6];
    int total_lines=0;
    string line,word;
    int lines;
    ifstream file;
    file.open(file_name);
    if (!file.is_open())
    {
        cout << "Can't open dump.atom file." << endl;
        exit(1); 
    }

//comfirm atom_number:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    for (int i=0;i<3;i++)
    {
        getline(file,line);
    }
    file >> atom_number;   //read the total atom number
    getline(file,line);
    getline(file,line);
    for (int i=0;i<3;i++)
    {
        getline(file,line);
        istringstream WORDS(line);
        for (int j=0;j<2;j++)
        {
            if (WORDS >> word)
            {
                if (j==0)
            {
                x_b[i*2+j]=stod(word);
            }
            else if (j==1)
            {
                x_b[i*2+j]=stod(word);
            }
            }
        }
    }
    a=x_b[1]-x_b[0];
    b=x_b[3]-x_b[2];
    c=x_b[5]-x_b[4];
    getline(file,line);
    istringstream WORDS(line);
    if (WORDS >> word)
    {
        if (word=="xy")
            lines=10;
        else
            lines=9;
    }
    file.clear();
    file.seekg(0,ios::beg);
//comfirm steps:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    while(getline(file,line))
    {
        ++total_lines;
    }
    step=total_lines/(lines+atom_number);
    file.clear();
    file.seekg(0,ios::beg);
//read xyz:!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    for (int i=0;i<step;i++)
    {
        for (int j=0;j<lines;j++)
        {
            getline(file,line);
        }

        for (int j=0;j<atom_number;j++)
        {
            getline(file,line);
            istringstream words(line);
            int indices,type;
            double x,y,z;
            for (int k=0;k<5;k++)
            {
                if (words >> word)
                {
                    if (k==0)
                    {
                        indices=stoi(word);
                    }
                    else if (k==1)
                    {
                        type=stoi(word)-1;
                    }
                    else if (k==2)
                    {
                        x=stod(word);
                    }
                    else if (k==3)
                    {
                        y=stod(word);
                    }
                    else if (k==4)
                    {
                        z=stod(word);
                    }
                }
            }
            R_old.push_back(atom_old(i,indices,type,x,y,z));
        }
        if (i%static_cast<int>(round(0.2*step))==0)
        {
            cout << i << " is read. " << endl;
        }
    }
    sort(R_old.begin(),R_old.end());
}
//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                         moving xyz                                   //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
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
                if (dx>=(0.5*a))
                {
                    d_R[(i+1)*atom_number+n].x+=a;
                }
                else if (dx<=(-0.5*a))
                {
                    d_R[(i+1)*atom_number+n].x-=a;
                }

                if (dy>=(0.5*b))
                {
                    d_R[(i+1)*atom_number+n].y+=b;
                }
                else if (dy<=(-0.5*b))
                {
                    d_R[(i+1)*atom_number+n].y-=b;
                }

                if (dz>=(0.5*c))
                {
                    d_R[(i+1)*atom_number+n].z+=c;
                }
                else if (dz<=(-0.5*c))
                {
                    d_R[(i+1)*atom_number+n].z-=c;
                }
            }
        }
    } 
}
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
//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                         MSD_WRITE                                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
void MSD_write(double *h_msd_ave)
{
    ofstream write_file;
    write_file.open("data/msd.dat");
    write_file.precision(9);
    write_file.setf(ios::fixed, ios::floatfield);
    if (!write_file.is_open())
    {
        cout << "Can't open dump.atom file." << endl;
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