msd_gpu: msd.cu
	nvcc -O3 -g -std=c++11 msd.cu -o msd_gpu \
	-gencode arch=compute_70,code=sm_70 \
	-gencode arch=compute_80,code=sm_80
