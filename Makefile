SRC_DIR = src

SRCS = $(SRC_DIR)/msd.cu $(SRC_DIR)/boundary.cu $(SRC_DIR)/calculation.cu $(SRC_DIR)/read_dump.cpp $(SRC_DIR)/write_dump.cpp

TARGET = msd_gpu

NVCC = nvcc
NVCC_FLAGS = -O3 -g -std=c++11
GENCODE_FLAGS = -gencode arch=compute_70,code=sm_70 -gencode arch=compute_80,code=sm_80

all: $(TARGET)

$(TARGET): $(SRCS)
	$(NVCC) $(NVCC_FLAGS) $(GENCODE_FLAGS) $^ -o $@

clean:
	rm -f $(TARGET)
