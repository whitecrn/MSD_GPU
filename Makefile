SRC_DIR = src

SRCS = $(SRC_DIR)/msd.cu $(SRC_DIR)/boundary.cu $(SRC_DIR)/calculation.cu $(SRC_DIR)/read_dump.cpp $(SRC_DIR)/write_msd.cpp $(SRC_DIR)/information.cu

TARGET = msd_gpu

CUDA_PATH = /usr/local/cuda
NVCC = $(CUDA_PATH)/bin/nvcc
NVCC_FLAGS = -O3 -g -std=c++11
GENCODE_FLAGS = -gencode arch=compute_70,code=sm_70 -gencode arch=compute_80,code=sm_80

all: $(TARGET)

$(TARGET): $(SRCS)
	$(NVCC) $(NVCC_FLAGS) $(GENCODE_FLAGS) $^ -o $@

clean:
	rm -f $(TARGET)
