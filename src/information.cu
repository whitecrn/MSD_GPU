#include "../include/information.h"

int show_information()
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
        return 1;
    }

    cout << "Device " << device_id << ": " << device_prop.name << endl;
    cout << "Compute capability: " << device_prop.major << "." << device_prop.minor << endl;
    return 0;
}