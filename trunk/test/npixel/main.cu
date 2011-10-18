#include <assert.h>
#include <stdio.h>

#define PIXEL_SIZE (512)
#define PIXEL_BLOCK	(16)
#define TILE_SIZE (16)
#define GRID_SIZE (PIXEL_SIZE/TILE_SIZE)

int total_r(unsigned int *results, unsigned int size){
  unsigned int i, total = 0;
  for(i = 0; i < size; ++i){
        total += results[i];
  }
  return total;
}

__global__ void glensing(const float pixelBlock, unsigned int* results) {
	
  	const int col = blockDim.x * blockIdx.x;
	const int row = pixelBlock * PIXEL_BLOCK + blockDim.y * blockIdx.y;
    
    atomicAdd(&results[row * PIXEL_SIZE + col], 1);
}

int main(){
	unsigned int *result = (unsigned int*)malloc(PIXEL_SIZE * PIXEL_SIZE * sizeof(unsigned int));
	unsigned int *result_d;
    cudaMalloc(&result_d, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
	cudaMemset(result_d, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));

	dim3 gdim(PIXEL_SIZE, PIXEL_BLOCK);
  	dim3 bdim(TILE_SIZE, TILE_SIZE);
  	int proc_row;
  
	for(proc_row=0; proc_row< PIXEL_SIZE/PIXEL_BLOCK; ++proc_row){
		glensing<<<gdim, bdim>>>(proc_row, result_d);
		cudaThreadSynchronize();  
  	}
  	cudaMemcpy(result, result_d, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	
	int total = total_r(result, PIXEL_SIZE * PIXEL_SIZE);
	printf("the number of total rays is %d\n", total);

	
	cudaFree(result_d);
	
	return 0;	
}
