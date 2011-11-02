#include "global.h"
#include "util.h"
#include "constants.h"

#include <assert.h>
#include <stdio.h>
#include <omp.h>

/*selecting one of the strategies on threads management*/
/*setting the number after Pixel_Per_Block as '1' will allocate pixel per block, if the number is zero then allocate pixel per thread*/
#define PIXEL_PER_BLOCK		(0)
#define	PIXEL_PER_THREAD	(!PIXEL_PER_BLOCK)

/*setting up the size of the image pixel*/
#define PIXEL_SIZE 	(512)

/*the number of the size of the block mattters for PIXEL_PER_BLOCK, the number here is suppsoed to be optimal*/
#define BLOCK_SIZE	(8)

/*the tile size and the grid size, matters when the threads are managed as PIXEL_PER_THREAD, the numbers here are supposed to be optimal*/
#define TILE_SIZE	(16)
#define GRID_SIZE 	(PIXEL_SIZE/TILE_SIZE)

/*put the number of GPUs available here*/
#define NUM_DEVICE	(1)

/*Pointers to the lens x,y co-ordinates and mass on the global memory*/
float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;

/*the struct that contains some key constant values that will be passed into GPU*/
typedef struct d_constants{
  unsigned int rpp;
  float kappa_c, gamma_, source_scale;
  float image_scale_x, image_scale_y;
  float increment_x, increment_y;
} vars;

/*setting up some key parameters for the */
void init_variables(d_constants *const_struct) {
  const_struct->rpp = 32;
  const_struct->kappa_c = kappa_c;
  const_struct->gamma_ = gamma_;
  const_struct->source_scale = source_scale;
  const_struct->image_scale_x = image_scale_x;
  const_struct->image_scale_y = image_scale_y;
  const_struct->increment_x = 0;
  const_struct->increment_y = 0;
}

int highest(unsigned int *results, unsigned int size) {
  unsigned int i, highest_count = 0;
  for(i = 0; i < size; ++i){
    if (results[i] > highest_count)
      highest_count = results[i];
  }
  return highest_count;
}


int main(int argc, char** argv){
  
  float increment_x, increment_y;
  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  d_constants *const_struct = (d_constants*)salloc(sizeof(d_constants));
  init_variables(const_struct);
  read_lenses(argv[1]);
  
  
  int num_devices;
  cudaGetDeviceCount (&num_devices);
  
  if(num_devices != NUM_DEVICE) error("Wrong configuration on the number of devices, please re-confirm the number of devices");
  
  fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
  increment_x = (image_scale_x * 2) / PIXEL_SIZE;
  increment_y = (image_scale_y * 2) / PIXEL_SIZE;
  
  const_struct->increment_x = increment_x;
  const_struct->increment_y = increment_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);

  unsigned int *results[NUM_DEVICE];
  int i;
  for(i=0; i<NUM_DEVICE; ++i)
  results[i] = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
  if (!results) error("calloc failed in allocating the result array");
  
    omp_set_num_threads(NUM_DEVICE);  // create as many CPU threads as there are CUDA devices
  #pragma omp parallel
  {
    	/* Pointers to the lens x,y co-ordinates and mass on the GPU device */
		float *d_lens_x;
		float *d_lens_y;
		float *d_lens_mass;
		unsigned int *d_results;
		int device_No = omp_get_thread_num();
		d_constants *d_const_struct;
		
		cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
	  	cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
	  	cudaMalloc(&d_lens_mass, sizeof(float) * nobjects);
	 	cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
	  	cudaMalloc(&d_const_struct, sizeof(d_constants));
	
	  	cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
	  	cudaMemcpy(d_const_struct, const_struct, sizeof(vars), cudaMemcpyHostToDevice);
	  	cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
	  	cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
	  	cudaMemcpy(d_lens_mass, lens_mass, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
	  	
	  	
	  	/*free the variables*/
	  	cudaFree(d_lens_x);	  
	  	cudaFree(d_lens_y);
	  	cudaFree(d_lens_mass);
	  	cudaFree(d_results);
		cudaFree(d_const_struct);    	
  }
  /*free variables*/
  free(lens_x);
  free(lens_y);
  free(lens_mass);
  free(results);
  free(const_struct);
   
  /*exit the program*/
  return 0;
}
