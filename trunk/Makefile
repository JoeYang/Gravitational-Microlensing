#! /usr/bin/env make

BINARY = grav_lensing

CFLAGS =  -g -G

CUDA_DIR = /usr/local/cuda
DEVELOPMENT_DIR = ./debug
PRODUCTION_DIR = ./build

NVCC = $(CUDA_DIR)/bin/nvcc
NVCC_CFLAGS = -I$(CUDA_DIR)/include --compiler-options -W
#NVCC_LDFLAGS = -L$(CUDA_DIR)/lib64

DIRS = $(DEVELOPMENT_DIR) $(PRODUCTION_DIR)
SOURCE_DIR = ./src
SOURCES = $(shell find $(SOURCE_DIR) -name "*.cu")
DEVELOPMENT_BINARY = $(DEVELOPMENT_DIR)/$(BINARY)
PRODUCTION_BINARY = $(PRODUCTION_DIR)/$(BINARY)


.PHONY: all clean depends production development

# phony targets
all: main

main: $(DIRS) $(DEVELOPMENT_BINARY) $(PRODUCTION_BINARY)
	
production: $(PRODUCTION_DIR) $(PRODUCTION_BINARY)	

development: $(DEVELOPMENT_DIR) $(DEVELOPMENT_BINARY)

clean:
	rm -rf $(DEVELOPMENT_DIR) $(DEVELOPMENT_BINARY)
	rm -rf $(PRODUCTION_DIR) $(PRODUCTION_BINARY)

depends:
	$(NVCC) $(NVCC_CFLAGS) -M $(SOURCES) > Makefile.depends
	
$(DEVELOPMENT_DIR)/%.o: %.cu
	$(NVCC) $(NVCC_CFLAGS) -c -o $@ $<

$(DEVELOPMENT_BINARY): $(SOURCES)
	$(NVCC) $(NVCC_LDFLAGS) -o $@ $^
	
$(DEVELOPMENT_DIR):
	mkdir -p $(DEVELOPMENT_DIR)

$(PRODUCTION_DIR)/%.o: %.cu
	$(NVCC) $(NVCC_CFLAGS) -c -o $@ $<

$(PRODUCTION_BINARY): $(SOURCES)
	$(NVCC) $(NVCC_LDFLAGS) -o $@ $^
	
$(PRODUCTION_DIR):
	mkdir -p $(PRODUCTION_DIR)
		
# dependencies
-include Makefile.depends

