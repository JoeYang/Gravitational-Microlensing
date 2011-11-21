extern "C" {
#include "global.h"
#include "util.h"
#include "constants.h"
}
#include <assert.h>
#include <stdio.h>
#include <omp.h>
#include <cuda.h>
//#include <curand.h>
//#include <curand_kernel.h>

/*selecting one of the strategies on threads management*/
/*setting the number after Pixel_Per_Block as '1' will allocate pixel per block, if the number is zero then allocate pixel per thread*/
#define PIXEL_PER_BLOCK		(0)
#define	PIXEL_PER_THREAD	(!PIXEL_PER_BLOCK)

/*put the number of GPUs available here*/
#define NUM_DEVICE	(2)

#define KERNEL_CALL_NUM	(4)

#if PIXEL_PER_BLOCK
	#define RPP	1
#else
	#define RPP (32)
#endif		

/*setting up the size of the image pixel*/
#define PIXEL_SIZE 	(1024)

/*the number of the size of the block mattters for PIXEL_PER_BLOCK, the number here is suppsoed to be optimal*/
#define BLOCK_SIZE	(8)

/*the tile size and the grid size, matters when the threads are managed as PIXEL_PER_THREAD, the numbers here are supposed to be optimal*/
#define TILE_SIZE	(16)
#define GRID_SIZE 	(PIXEL_SIZE/TILE_SIZE)

/*Pointers to the lens x,y co-ordinates and mass on the global memory*/
float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;

/*setting up some key parameters for the */
void init_variables(d_constants *const_struct) {
  const_struct->rpp = RPP;
  const_struct->kappa_c = kappa_c;
  const_struct->gamma_ = gamma_;
  const_struct->source_scale = source_scale;
  const_struct->image_scale_x = image_scale_x;
  const_struct->image_scale_y = image_scale_y;
  const_struct->increment_x = 0;
  const_struct->increment_y = 0;
}

__global__ void glensing_rpt(const float *lens_x, const float *lens_y, const float *lens_mass, const size_t nobjects, unsigned int* results, const d_constants* v) {
  const unsigned int row = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int col = blockIdx.y*blockDim.y + threadIdx.y;
  
  const float initial_x = (-v->image_scale_x) + row*v->increment_x;
  const float initial_y = (-v->image_scale_y) + col*v->increment_y;

  const unsigned int uniform_box = sqrtf((float)v->rpp);
  const float source_scale = v->source_scale;
 
  float start_x, start_y, dx, dy;
  unsigned int it, noise_x, noise_y;
  size_t k;
  
  for(it = 0; it < v->rpp; ++it) {
    noise_x = it % uniform_box;
    noise_y = it / uniform_box;
    start_x = initial_x + noise_x * v->increment_x / uniform_box;
    start_y = initial_y + noise_y * v->increment_y / uniform_box;

    dx = (1-v->gamma_)*start_x - v->kappa_c*start_x;
    dy = (1+v->gamma_)*start_y - v->kappa_c*start_y;
	
	for(k = 0; k < nobjects; ++k) {
      float dist = pow(start_x - lens_x[k], 2) + pow(start_y - lens_y[k], 2);
      dx -= lens_mass[k] * (start_x - lens_x[k]) / dist;
      dy -= lens_mass[k] * (start_y - lens_y[k]) / dist;
    }
 
    if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
        (dy >= -source_scale/2) && (dy <= source_scale/2)) {
      int px = (dx + source_scale/2) / (source_scale/PIXEL_SIZE);
      int py = PIXEL_SIZE - (dy + source_scale/2) / (source_scale/PIXEL_SIZE);
      atomicAdd(&results[py * PIXEL_SIZE + px], 1);
    }
  }
}

__global__ void glensing_rpb(const float *lens_x, float *lens_y, float *lens_mass, const size_t nobjects, unsigned int* results, const d_constants* v) {
  	const float col = blockIdx.x;
	const float row = blockIdx.y;
	const float bx = threadIdx.x;
	const float by = threadIdx.y;	
  	
  	//Position of each light ray inside each Block	
	__device__ __shared__ float unit_x;
	__device__ __shared__ float unit_y;
	__device__ __shared__ float source_scale;  
 	__device__ __shared__ float base_x;
  	__device__ __shared__ float base_y; 
  	
  	if((bx+by) == 0){
	  	base_x = (-v->image_scale_x) + row*v->increment_x;
	  	base_y = (-v->image_scale_y) + col*v->increment_y;
	  	unit_x = v->increment_x/TILE_SIZE;
	  	unit_y = v->increment_y/TILE_SIZE;	
		source_scale = v->source_scale;
  	}
  	__syncthreads();
  	
	float start_x, start_y, dx, dy;
	size_t k;
	float dist;

	start_x = base_x + (bx) * unit_x;
	start_y = base_y + (by) * unit_y; 
    dx = (1-v->gamma_)*start_x - v->kappa_c*start_x;
    dy = (1+v->gamma_)*start_y - v->kappa_c*start_y;

    for(k = 0; k < nobjects; ++k) {
      dist = pow(start_x - lens_x[k], 2) + pow(start_y - lens_y[k], 2);
      dx -= lens_mass[k] * (start_x - lens_x[k]) / dist;
      dy -= lens_mass[k] * (start_y - lens_y[k]) / dist;
    }
    	
    if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
        (dy >= -source_scale/2) && (dy <= source_scale/2)) {
     	 int px = (dx + source_scale/2) / (source_scale/PIXEL_SIZE);
     	 int py = PIXEL_SIZE - (dy + source_scale/2) / (source_scale/PIXEL_SIZE);    	 
    	 atomicAdd(&results[py * PIXEL_SIZE + px], 1);
    }
}


int main(int argc, char** argv){
  
  float increment_x, increment_y;
  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  d_constants *const_struct = (d_constants*)salloc(sizeof(d_constants));
  init_variables(const_struct);
  read_lenses(argv[1]);
  
  fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
  increment_x = (image_scale_x * 2) / PIXEL_SIZE;
  increment_y = (image_scale_y * 2) / PIXEL_SIZE;
  const_struct->increment_x = increment_x;
  const_struct->increment_y = increment_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);
  
  int num_devices;
  cudaGetDeviceCount (&num_devices); 
  if(num_devices != NUM_DEVICE) error("Wrong configuration on the number of devices, please re-confirm the number of devices");
   
	printf("---------------------------\n");
	printf("%d GPGPU device found:\n", num_devices);
	for(int i = 0; i < num_devices; i++)
	{
		cudaDeviceProp dprop;
		cudaGetDeviceProperties(&dprop, i);
				printf("   %d: %s\n", i, dprop.name);
	}
  printf("---------------------------\n");

  unsigned int *results[KERNEL_CALL_NUM];
  int i;
  for(i=0; i<KERNEL_CALL_NUM; ++i){
  	results[i] = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
  	if (!results[i]) error("calloc failed in allocating the result array");
  }

	omp_set_num_threads(KERNEL_CALL_NUM);	
#pragma omp parallel
	{	
		unsigned int num_cpu_threads = omp_get_num_threads();
		unsigned int CPU_thread_No = omp_get_thread_num();
		cudaSetDevice(CPU_thread_No%NUM_DEVICE);
		int device_No = -1;
		cudaGetDevice(&device_No);
		printf("This is thread %d (of %d), using %d device\n", CPU_thread_No, num_cpu_threads, device_No);
    	
		
		/* Pointers to the lens x,y co-ordinates and mass on the GPU device */
		float *d_lens_x;
		float *d_lens_y;
		float *d_lens_mass;
		unsigned int *d_results;
		d_constants *d_const_struct;
		
		cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
	  	cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
	  	cudaMalloc(&d_lens_mass, sizeof(float) * nobjects);
	 	cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
	  	cudaMalloc(&d_const_struct, sizeof(d_constants));
	
	  	cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
	  	cudaMemcpy(d_const_struct, const_struct, sizeof(d_constants), cudaMemcpyHostToDevice);
	  	cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
	  	cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
	  	cudaMemcpy(d_lens_mass, lens_mass, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
	  	
	  	if(PIXEL_PER_BLOCK){	  		
	  	  // Perform gravitational microlensing
		  dim3 bdim(BLOCK_SIZE, BLOCK_SIZE);
		  dim3 gdim(PIXEL_SIZE, PIXEL_SIZE);
		  glensing_rpb<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_const_struct);
		  cudaMemcpy(results[CPU_thread_No], d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost);
	  	}
	  	else if(PIXEL_PER_THREAD){	  		
		  // Perform gravitational microlensing
		  dim3 bdim(TILE_SIZE, TILE_SIZE);
		  dim3 gdim(GRID_SIZE, GRID_SIZE);  
		  glensing_rpt<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_const_struct);
  	  	  cudaMemcpy(results[CPU_thread_No], d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost);

	  	}
	  	
	  	/*free the variables*/
	  	cudaFree(d_lens_x);	  
	  	cudaFree(d_lens_y);
	  	cudaFree(d_lens_mass);
	  	cudaFree(d_results);
		cudaFree(d_const_struct);
	}
  
  unsigned int *final_result = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
	
	int r_c=0, t;
	for(; r_c < PIXEL_SIZE*PIXEL_SIZE; ++r_c){
		for(t=0; t<KERNEL_CALL_NUM; ++t)
			final_result[r_c] += results[t][r_c];
	}
	
	int total = total_r(final_result, PIXEL_SIZE * PIXEL_SIZE);
	printf("The total num of rays is %d\n", total);

   int highest_c = highest(final_result, PIXEL_SIZE * PIXEL_SIZE);
   write_pgm(final_result, PIXEL_SIZE, PIXEL_SIZE, highest_c);
  
  /*free variables*/
  free(lens_x);
  free(lens_y);
  free(lens_mass);
  free(const_struct);
  free(final_result);
  for(i=0; i<KERNEL_CALL_NUM; ++i)
  	free(results[i]);
   
  cudaThreadExit();
  /*exit the program*/
  return 0;
}
