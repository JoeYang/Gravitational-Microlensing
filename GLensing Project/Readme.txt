This is the project for course INFO5011.

Members:
	Xiangrui Yang, Katrina Nally, Andrew Stevens, Stephen Merity

In this directory there are 6 folders, containing the source code for six different implementations.

CPU - the CPU code with sequential calculation.
GPU - Naive parallism implementation
OpenMP 1.0 - the first edition of openMP and CUDA implementation
OpenMP 2.0 - the second edition of openMP and CUDA implementation
cudaTreeSrc - tree optimization with CUDA
OpenMP Tree - openMP with tree optimisation

When trying to execute the program, please make sure that the pixel size, block size, number of omp threads are created properly. Please try to make a rough approximation about the memory usage and make sure that the graphical memory is enough for all the processes.
