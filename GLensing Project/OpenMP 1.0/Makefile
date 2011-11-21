NVCCFLAGS	:= -O3 -g -lm
NVCC		:= /usr/local/cuda/bin/nvcc
LD_LIBRARY_PATH	:= -L/usr/local/cuda/lib64 

all: glensing

glensing: main.cu constants.c util.c 
	$(NVCC) $(NVCCFLAGS) $(LD_LIBRARY_PATH) -arch compute_11 -Xcompiler -fopenmp -o $@ $^
	
clean:
	rm -rf *.o
	rm -rf glensing
