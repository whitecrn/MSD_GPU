msd_gpu: msd.cu boundary.cu calculation.cu read_dump.cpp write_dump.cpp
	nvcc -O3 -g -std=c++11 msd.cu boundary.cu calculation.cu read_dump.cpp write_dump.cpp -o msd_gpu \
-gencode arch=compute_70,code=sm_70 \
-gencode arch=compute_80,code=sm_80 \
