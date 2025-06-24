#include "../include/msd.h"
#include "../include/type.h"

double potim=0.1;
int atom_id=0;
int step=0;
string file_name="dump.atom";
bool h_flag=false;
bool Cartesian=true;

int main(int argc,char *argv[])
{   
    int opt;
    while ((opt=getopt(argc,argv,"i:p:f:t:h"))!=-1)
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
            case 't':
                Cartesian=atoi(optarg);
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