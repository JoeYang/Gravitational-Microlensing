# GPGPU
BINARY = gpgpu_glens

PRODUCTION_CFLAGS = -O3
PRODUCTION_DIR = build
DEVELOPMENT_CFLAGS = -g --ptxas-options=-v
DEVELOPMENT_DIR = debug

# Set up the directories for CUDA on the GPU machines
CUDA_DIR = /usr/local/cuda

NVCC = $(CUDA_DIR)/bin/nvcc
NVCC_CFLAGS = -I$(CUDA_DIR)/include --compiler-options -W -Xopencc -OPT:unroll_size=200000 -G
NVCC_LDFLAGS = -L$(CUDA_DIR)/open64 -lcurand  -Xcompiler -fopenmp  -arch compute_11

SRC_DIR = ./src
GPU_DIRS = $(PRODUCTION_DIR) $(DEVELOPMENT_DIR)
GPU_SRC = $(SRC_DIR)/main.cu $(SRC_DIR)/util.h $(SRC_DIR)/constants.h $(SRC_DIR)/global.h $(SRC_DIR)/constants.c $(SRC_DIR)/util.c
GPU_BIN = $(PRODUCTION_DIR)/$(BINARY)
GPU_DEV_BIN = $(DEVELOPMENT_DIR)/$(BINARY)

.PHONY: all clean depends production development

# phony targets
all: $(GPU_DIRS) $(GPU_BIN) $(GPU_DEV_BIN)

production: $(PRODUCTION_DIR) $(PRODUCTION_BINARY)

development: $(DEVELOPMENT_DIR) $(DEVELOPMENT_BINARY)

clean:
	-rm -f $(GPU_BIN) $(GPU_OBJ)
	-rm -f $(GPU_DEV_BIN) $(GPU_DEV_OBJ)

depends:
	$(NVCC) $(NVCC_CFLAGS) -M $(GPU_SRC) > Makefile.depends

# real targets
$(PRODUCTION_DIR)/%.o: %.cu
	$(NVCC) $(PRODUCTION_CFLAGS) -c -o $@ $<
$(PRODUCTION_DIR)/%.o: %.c
	$(NVCC) $(PRODUCTION_CFLAGS) -c -o $@ $<

$(DEVELOPMENT_DIR)/%.o: %.cu
	$(NVCC) $(DEVELOPMENT_CFLAGS) -c -o $@ $<
$(DEVELOPMENT_DIR)/%.o: %.c
	$(NVCC) $(DEVELOPMENT_CFLAGS) -c -o $@ $<

$(GPU_BIN): $(GPU_SRC)
	$(NVCC) $(NVCC_LDFLAGS) -o $@ $^
$(GPU_DEV_BIN): $(GPU_SRC)
	$(NVCC) $(NVCC_CFLAGS) $(NVCC_LDFLAGS) -o $@ $^

$(DEVELOPMENT_DIR):
	mkdir -p $(DEVELOPMENT_DIR)
$(PRODUCTION_DIR):
	mkdir -p $(PRODUCTION_DIR)


