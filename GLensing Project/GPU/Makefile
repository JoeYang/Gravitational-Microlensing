# GPGPU
BINARY = gpgpu_glens

PRODUCTION_CFLAGS = -arch compute_11
PRODUCTION_DIR = build
DEVELOPMENT_CFLAGS = -arch compute_11 -g -G -deviceemu
DEVELOPMENT_DIR = debug

# Set up the directories for CUDA on the GPU machines
ifeq ($(shell hostname | grep -q '^gpu' && echo yes || echo no), yes)
	CUDA_DIR = /usr/local/cuda
else
	CUDA_DIR = /opt/cuda
endif

NVCC = $(CUDA_DIR)/bin/nvcc
NVCC_CFLAGS = -I$(CUDA_DIR)/include --compiler-options -W
NVCC_LDFLAGS = -L$(CUDA_DIR)/lib64

GPU_DIRS = $(PRODUCTION_DIR) $(DEVELOPMENT_DIR)
GPU_SRC = main.cu util.c util.h constants.c constants.h
GPU_OBJ = $(PRODUCTION_DIR)/main.o $(PRODUCTION_DIR)/util.o $(PRODUCTION_DIR)/constants.o
GPU_DEV_OBJ = $(DEVELOPMENT_DIR)/main.o $(DEVELOPMENT_DIR)/util.o $(DEVELOPMENT_DIR)/constants.o
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
	$(NVCC) $(NVCC_CFLAGS) $(PRODUCTION_CFLAGS) -c -o $@ $<
$(PRODUCTION_DIR)/%.o: %.c
	$(NVCC) $(NVCC_CFLAGS) $(PRODUCTION_CFLAGS) -c -o $@ $<

$(DEVELOPMENT_DIR)/%.o: %.cu
	$(NVCC) $(NVCC_CFLAGS) $(DEVELOPMENT_CFLAGS) -c -o $@ $<
$(DEVELOPMENT_DIR)/%.o: %.c
	$(NVCC) $(NVCC_CFLAGS) $(DEVELOPMENT_CFLAGS) -c -o $@ $<

$(GPU_BIN): $(GPU_OBJ)
	$(NVCC) $(NVCC_LDFLAGS) -o $@ $^
$(GPU_DEV_BIN): $(GPU_DEV_OBJ)
	$(NVCC) $(NVCC_LDFLAGS) -o $@ $^

$(DEVELOPMENT_DIR):
	mkdir -p $(DEVELOPMENT_DIR)
$(PRODUCTION_DIR):
	mkdir -p $(PRODUCTION_DIR)

# dependencies
-include Makefile.depends
