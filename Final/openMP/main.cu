extern "C"{
#include "global.h"
#include "util.h"
#include "constants.h"
#include "tree_struct.h"
}
#include <assert.h>
#include <stdio.h>
#include <omp.h>
#include <cuda.h>
#include "curand.h"
#include "curand_kernel.h"

#define PIXEL_SIZE 	(1024)
#define pixel_x (PIXEL_SIZE)
#define pixel_y (PIXEL_SIZE)
#define QUAD_IDX (4)

/*put the number of GPUs available here*/
#define NUM_DEVICE	(2)

#define KERNEL_CALL_PER_DEVICE	(4)
#define KERNEL_CALL_NUM	(KERNEL_CALL_PER_DEVICE * NUM_DEVICE)

//number of rays per pixel
#define RPP (512)
#define ITERATION	(RPP/KERNEL_CALL_NUM)
/*setting up the size of the image pixel*/

/*the tile size and the grid size, matters when the threads are managed as PIXEL_PER_THREAD, the numbers here are supposed to be optimal*/
#define TILE_SIZE	(16)
#define GRID_SIZE 	(PIXEL_SIZE/TILE_SIZE)

float4 *lenses, *d_lenses;
/*Pointers to the lens x,y co-ordinates and mass on the global memory*/
float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;
cell *root; 	/* global pointer to root cell so we can get to it from anywhere */
const float delta = 0;		/*variable to determine accuracy*/
static int lens_index = 0;

/*setting up some key parameters for the */
void init_variables(d_constants *const_struct) {
  const_struct->rpp = 1;
  const_struct->kappa_c = kappa_c;
  const_struct->gamma_ = gamma_;
  const_struct->source_scale = source_scale;
  const_struct->image_scale_x = image_scale_x;
  const_struct->image_scale_y = image_scale_y;
  const_struct->increment_x = 0;
  const_struct->increment_y = 0;
}

/*
* Method to get the required lenses for the lensing calculation
*/
void get_lenses(float4 ** lenses, cell ** tree, float delta, float ray_x, float ray_y){
  int j, i, has_cells=0, has_lenses=0;

  /* finding the values of the current cell to decide if it is to be included in the calculation*/
  float cell_x = (float)(*tree)->center_mass_x/(float)(*tree)->total_mass;
  float cell_y = (float)(*tree)->center_mass_y/(float)(*tree)->total_mass;
  float dist = sqrt(pow(ray_x - cell_x,2) + pow(ray_y - cell_y, 2));
  float ratio = ((*tree)->width)/dist;

  // checking if the current cell has any cells or lenses
  for(j=0; j<QUAD_IDX; j++){
    if((*tree)->lenses[j]) has_lenses++;
    if((*tree)->desc[j]) has_cells++;
  }

  if(delta == 0.0 && has_lenses>0){
    //lens_index += has_lenses;
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->lenses[i]){
        (*lenses)[lens_index].x = (*tree)->lenses[i]->x;
        (*lenses)[lens_index].y = (*tree)->lenses[i]->y;
        (*lenses)[lens_index].w = (*tree)->lenses[i]->m;
        lens_index++;
      }
      if((*tree)->desc[i]) get_lenses(lenses, &((*tree)->desc[i]), delta, ray_x, ray_y);
    }
  }
  else if(has_lenses==1 && has_cells==0){
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->lenses[i]){
        (*lenses)[lens_index].x = (*tree)->lenses[i]->x;
        (*lenses)[lens_index].y = (*tree)->lenses[i]->y;
        (*lenses)[lens_index].w = (*tree)->lenses[i]->m;
      }
    }
    lens_index++;
    return;
  }
  else if(ratio<delta){
    (*lenses)[lens_index].x = cell_x;
    (*lenses)[lens_index].y = cell_y;
    (*lenses)[lens_index].w =(*tree)->total_mass;
    lens_index++;
    return;
  }
  else{
    int i;
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->desc[i]) get_lenses(lenses, &((*tree)->desc[i]), delta, ray_x, ray_y);
    }
  }
}

__global__ void curand_setup(curandState* globalState, unsigned long seed){
	const unsigned int row = blockIdx.x;
	const unsigned int col = blockIdx.y;
	const unsigned int id = row*PIXEL_SIZE + col;
	curand_init(seed, id, 0, &globalState[id]);
}

__global__ void glensing(const float4 *lenses,  const size_t nobjects, unsigned int* results, const d_constants* v, curandState* globalState){
	const unsigned int row = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int col = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int id = row*PIXEL_SIZE + col;
	const unsigned int lens_idx = blockIdx.y*blockDim.y + blockIdx.x*nobjects;
	const float initial_x = (-v->image_scale_x) + row*v->increment_x;
	const float initial_y = (-v->image_scale_y) + col*v->increment_y;
	const float increment_x = v->increment_x;
	const float increment_y = v->increment_y;
	const float source_scale = v->source_scale;
	__device__ __shared__ curandState localState;
	
	
	float start_x, start_y, dx, dy;
	size_t k;
	localState = globalState[id];
	start_x = initial_x + curand_uniform(&localState) * increment_x;
	globalState[id] = localState;
	start_y = initial_y + curand_uniform(&localState) * increment_y;
	globalState[id] = localState;
	
	dx = (1-v->gamma_)*start_x - v->kappa_c*start_x;
	dy = (1+v->gamma_)*start_y - v->kappa_c*start_y;

	for(k = 0; k < nobjects; ++k) {
		float dist = pow(start_x - lenses[lens_idx + k].x, 2) + pow(start_y - lenses[lens_idx + k].y, 2);
		dx -= lenses[lens_idx + k].w * (start_x - lenses[lens_idx + k].x) / dist;
		dy -= lenses[lens_idx + k].w * (start_y - lenses[lens_idx + k].y) / dist;
	}

	if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
		(dy >= -source_scale/2) && (dy <= source_scale/2)) {
		int px = (dx + source_scale/2) / (source_scale/PIXEL_SIZE);
		int py = PIXEL_SIZE - (dy + source_scale/2) / (source_scale/PIXEL_SIZE);
		atomicAdd(&results[py * PIXEL_SIZE + px], 1);
	}
}

int main(int argc, char** argv){
	root = NULL;
	lens *temp_lens;
	
	int i, j, m, array_size;;
	float increment_x, increment_y;
	// Load relevant settings and data
	if (argc < 2) error("Requires argument with lens positions and optional mass");
	setup_constants();
	d_constants *const_struct = (d_constants*)salloc(sizeof(d_constants));
	init_variables(const_struct);
	read_lenses(argv[1]);
	make_root(&root);
	
	/* creating the tree*/
	for(i=0; i < (int)nobjects; i++){
		//temp = (cell *)salloc(sizeof(cell));                /*allocating memory to temp cell*/
		temp_lens = (lens *)salloc(sizeof(lens));                /*allocating memory to temp lens*/
		temp_lens->x = lens_x[i];                          /*assigning x-coord to cell_lens x-coord*/
		temp_lens->y = lens_y[i];                          /*assigning y-coord to cell_lens y-coord*/
		temp_lens->m = lens_mass[i];                          /*assigning mass*/
		make_tree(&root, temp_lens);
		free(temp_lens);
	}
	
	fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
	increment_x = (image_scale_x * 2) / (PIXEL_SIZE);
	increment_y = (image_scale_y * 2) / (PIXEL_SIZE);
	const_struct->increment_x = increment_x;
	const_struct->increment_y = increment_y;
	fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);
	
	int *num_lenses = (int *)salloc(sizeof(int));
	if((PIXEL_SIZE % TILE_SIZE) == 0){
		array_size = (PIXEL_SIZE/TILE_SIZE)*(PIXEL_SIZE/TILE_SIZE);
	}
	else{
		error("PIXEL_SIZE is not a multiple of TILE_SIZE");
	}
	
	*num_lenses = 0;
	float temp_x = -image_scale_x + increment_x;
	float temp_y = image_scale_y - increment_y;
	get_lens_count(&root, delta, temp_x, temp_y, num_lenses);
	lenses = (float4 *)salloc(sizeof(float4)*array_size*(*num_lenses));
	
  int count = 0;
/* performing all walks here on the CPU then transferring these to the GPU*/
	for(j=0; j<PIXEL_SIZE/TILE_SIZE; j++){
		for(i=0; i< PIXEL_SIZE/TILE_SIZE; i++){
		  *num_lenses = 0;
		  lens_index = 0;
		  /* we use these as we need the center coordinates of each tile to determine the lenses to be used for that tile*/
		  float tile_x = -image_scale_x + increment_x;
		  float tile_y = image_scale_y - increment_y;
		  get_lens_count(&root, delta, tile_x, tile_y, num_lenses);
		  float4 *inc_lenses = (float4 *)salloc(sizeof(float4)*(*num_lenses));
		  get_lenses(&inc_lenses, &root, delta, tile_x, tile_y);

		  for(m = 0; m<(*num_lenses); m++){
			int idx = j*(PIXEL_SIZE/TILE_SIZE)*(*num_lenses) +i*(*num_lenses) +m;
			lenses[idx] = inc_lenses[m];
		  }
		  free(inc_lenses);
		  count++;
		}
	}
	
	
	int num_devices;
	cudaGetDeviceCount (&num_devices); 
	if(num_devices != NUM_DEVICE) error("Wrong configuration on the number of devices, please re-confirm the number of devices");

	printf("---------------------------\n");
	printf("%d GPGPU device found:\n", num_devices);
	for(i = 0; i < num_devices; i++)
	{
		cudaDeviceProp dprop;
		cudaGetDeviceProperties(&dprop, i);
				printf("   %d: %s\n", i, dprop.name);
	}
	printf("---------------------------\n");

	unsigned int *results[KERNEL_CALL_NUM];

	for(i=0; i<KERNEL_CALL_NUM; ++i){
		results[i] = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
		if (!results[i]) error("calloc failed in allocating the result array");
	}
	
	unsigned int *final_result = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
	omp_set_num_threads(KERNEL_CALL_NUM);
	#pragma omp parallel 
	{
		int device_No = omp_get_thread_num()%NUM_DEVICE;
		unsigned int thread_No = omp_get_thread_num();
		cudaSetDevice(device_No);
		/* Pointers to the lens x,y co-ordinates and mass on the GPU device */

		d_constants *d_const_struct;
		curandState* globalState;
		
		cudaMalloc(&d_const_struct, sizeof(d_constants));
		cudaMalloc(&globalState, sizeof(curandState)*PIXEL_SIZE*PIXEL_SIZE/(TILE_SIZE*TILE_SIZE));
		
		dim3 bdim(TILE_SIZE, TILE_SIZE);
		dim3 gdim(GRID_SIZE, GRID_SIZE); 
		curand_setup<<<gdim, bdim>>>(globalState, time(NULL));
		cudaMemcpy(d_const_struct, const_struct, sizeof(d_constants), cudaMemcpyHostToDevice);
	
		cudaMalloc(&d_lenses, sizeof(float4)*(*num_lenses)*(PIXEL_SIZE/TILE_SIZE)*(PIXEL_SIZE/TILE_SIZE));
		cudaMemset(d_lenses, 0, sizeof(float4)*(*num_lenses)*(PIXEL_SIZE/TILE_SIZE)*(PIXEL_SIZE/TILE_SIZE));
		cudaMemcpy(d_lenses, lenses, sizeof(float4)*(*num_lenses)*(PIXEL_SIZE/TILE_SIZE)*(PIXEL_SIZE/TILE_SIZE), cudaMemcpyHostToDevice);
		
		int iteration;
		for(iteration = 0; iteration < ITERATION; ++iteration){		
			unsigned int *d_results;
			cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
			cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
			unsigned int *temp_results = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
			printf("Device: %d Thread: %d %4d/%4d Start\n", device_No, thread_No, iteration, ITERATION);
			// Perform gravitational microlensing
			glensing<<<gdim, bdim>>>(d_lenses, (*num_lenses), d_results, d_const_struct, globalState);
			//glensing<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_const_struct, globalState);
			cudaMemcpy(temp_results, d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost);  
			for(int k = 0; k<PIXEL_SIZE*PIXEL_SIZE; ++k){
				results[thread_No][k] += temp_results[k];
			}
			int total_x = total_r(temp_results, PIXEL_SIZE * PIXEL_SIZE);
			printf("Device: %d Thread: %d %4d/%4d Complete, total rays %ld\n", device_No, thread_No, iteration, ITERATION, total_x);
			free(temp_results);
			cudaFree(d_results);
		}		
		
		/*free the variables*/
		
		cudaFree(d_const_struct);
		cudaFree(globalState);
		cudaFree(d_lenses);
	}
	int r_c=0, t;
	for(r_c = 0; r_c < PIXEL_SIZE*PIXEL_SIZE; ++r_c){
		for(t=0; t<KERNEL_CALL_NUM; ++t)
			final_result[r_c] += results[t][r_c];
	}
	int total = total_r(final_result, PIXEL_SIZE * PIXEL_SIZE);
	printf("The total num of rays is %ld \n", total);

	int highest_c = highest(final_result, PIXEL_SIZE * PIXEL_SIZE);
	write_pgm(final_result, PIXEL_SIZE, PIXEL_SIZE, highest_c);
  
	/*free variables*/
	free(lens_x);
	free(lens_y);
	free(lens_mass);
	free(const_struct);
	free(final_result);
	free(lenses);
	free(num_lenses);
	free_tree(&root);
	for(i=0; i<KERNEL_CALL_NUM; ++i)
		free(results[i]);

	cudaThreadExit();
	/*exit the program*/
	return 0;
}
