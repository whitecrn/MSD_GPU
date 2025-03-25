#pragma once
#include <cstdio>
#include <iostream>
using namespace std;

#define CHECK(call) \
do \
{ \
    const cudaError_t error_code=call; \
    if (error_code!=cudaSuccess) \
    { \
        cout << "CUDA Error:\n" << endl; \
        cout << "  File:  " << __FILE__ << '\n' << endl; \
        cout << "  Line:  " << __LINE__ << '\n' << endl; \
        cout << "  Error code:  " << error_code << '\n' << endl; \
        cout << "  Error text:  " << cudaGetErrorString(error_code) << '\n' << endl; \
        exit(1); \
    } \
} while(0)
