// All C includes must be wrapped in extern "C"
extern "C" {
#include "global.h"
#include "util.h"
#include "constants.h"
}
#include <assert.h>
#include <stdio.h>

#define PIXEL_SIZE	128
#define TILE_SIZE 	4

float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;

/* Pointers to the lens x,y co-ordinates on the GPU device */
float *d_lens_x;
float *d_lens_y;
float *d_lens_mass;

typedef struct vars {
  unsigned int rpp;
  float kappa_c, gamma_, source_scale;
  float image_scale_x, image_scale_y;
  float increment_x, increment_y;
} vars;

void init_var(vars *var) {
  var->rpp = 1;
  var->kappa_c = kappa_c;
  var->gamma_ = gamma_;
  var->source_scale = source_scale;
  var->image_scale_x = image_scale_x;
  var->image_scale_y = image_scale_y;
  var->increment_x = 0;
  var->increment_y = 0;
}

int highest(unsigned int *results, unsigned int size) {
  unsigned int i, highest_count = 0;
  for(i = 0; i < size; ++i){
    if (results[i] > highest_count)
      highest_count = results[i];
  }
  return highest_count;
}

int total_r(unsigned int *results, unsigned int size){
  unsigned int i, total = 0;
  for(i = 0; i < size; ++i){
        total += results[i];
  }
  return total;
}

__global__ void glensing(const float *lens_x, float *lens_y, float *lens_mass, const size_t nobjects, unsigned int* results, const vars* v) {
	
  	const float col = blockIdx.x;
	const float row = blockIdx.y;
	
	const float bx = threadIdx.x;
	const float by = threadIdx.y; 

	__device__ __shared__ float base_x;
  	__device__ __shared__ float base_y;
  	
  	//Position of each light ray inside each Block	
	__device__ __shared__ float unit_x;
	__device__ __shared__ float unit_y;
	__device__ __shared__ float source_scale;  
 	
  	base_x = (-v->image_scale_x) + row*v->increment_x;
  	base_y = (-v->image_scale_y) + col*v->increment_y;
  	unit_x = v->increment_x/TILE_SIZE;
  	unit_y = v->increment_y/TILE_SIZE;	
	source_scale = v->source_scale;
	
	
	float start_x, start_y, dx, dy;
	size_t k;
	float dist;
   
    start_x = base_x + (bx ) * unit_x;
    start_y = base_y + (by ) * unit_y; 

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
	else
		atomicAdd(&results[0],1);
}

int main(int argc, char** argv) {  
	cudaSetDevice(1);
  float increment_x, increment_y;
  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  vars *variables = (vars *)salloc(sizeof(vars));
  init_var(variables);
  read_lenses(argv[1]);

  fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
  increment_x = (image_scale_x * 2) / PIXEL_SIZE;
  increment_y = (image_scale_y * 2) / PIXEL_SIZE;
  variables->increment_x = increment_x;
  variables->increment_y = increment_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);

  unsigned int *results = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
  unsigned int *d_results;
  if (!results) error("calloc failed in allocating the result array");

  // Setting up CUDA global memory
  vars *d_variables;
  cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
  cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
  cudaMalloc(&d_lens_mass, sizeof(float) * nobjects);
  
  cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
  cudaMalloc(&d_variables, sizeof(vars));

  cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
  cudaMemcpy(d_variables, variables, sizeof(vars), cudaMemcpyHostToDevice);
  cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lens_mass, lens_mass, sizeof(float) * nobjects, cudaMemcpyHostToDevice);

  // Perform gravitational microlensing
  dim3 bdim(TILE_SIZE, TILE_SIZE);
  dim3 gdim(PIXEL_SIZE, PIXEL_SIZE);
  
  glensing<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_variables);

   
  cudaMemcpy(results, d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost);

  int highest_c = highest(results, PIXEL_SIZE * PIXEL_SIZE);
  int total = total_r(results, PIXEL_SIZE * PIXEL_SIZE);
  write_pgm(results, PIXEL_SIZE, PIXEL_SIZE, highest_c);
  printf("the number of total rays is %d\n", total);
  
  
  // Free the memory allocated during processing
  // GPU
  cudaFree(d_results);
  cudaFree(d_variables);
  cudaFree(d_lens_x);
  cudaFree(d_lens_y);
  cudaFree(d_lens_mass);
  // CPU
  free(lens_x);
  free(lens_y);
  free(lens_mass);
  free(results);
  return 0;
}
