msd_gpu: msd.cu
	nvcc -O3 -g msd.cu -o msd_gpu \
	-gencode arch=compute_70,code=sm_70 \
	-gencode arch=compute_80,code=sm_80