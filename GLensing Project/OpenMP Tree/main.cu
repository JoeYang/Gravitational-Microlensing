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

#define PIXEL_SIZE 	(512)
#define QUAD_IDX (4)
/*put the number of GPUs available here*/
#define NUM_DEVICE	(2)
#define KERNEL_CALL_PER_DEVICE	(1)
#define KERNEL_CALL_NUM	(KERNEL_CALL_PER_DEVICE * NUM_DEVICE)
#define IF_TREE_OPTIMIZATION	1
#define TILE_SIZE	(16)
#define GRID_SIZE 	(PIXEL_SIZE/TILE_SIZE)

#define SEED_MODE	1024


//number of rays per pixel
#define RPP (16)
#define	RPP_Kernel	(1)
#define ITERATION	(RPP/(KERNEL_CALL_NUM*RPP_Kernel))

float *lenses_x, *lenses_y, *lenses_m;
float *d_lens_x, *d_lens_y, *d_lens_mass;
/*Pointers to the lens x,y co-ordinates and mass on the global memory*/
float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;
cell *root; 	/* global pointer to root cell so we can get to it from anywhere */
const float delta = 0.4;		/*variable to determine accuracy*/
static int lens_index = 0;

/*setting up some key parameters for the */
void init_variables(d_constants *const_struct) {
	const_struct->rpp = RPP_Kernel;
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
void get_lenses(float ** l_x, float ** l_y, float ** l_m, cell ** tree, float delta, float ray_x, float ray_y){
  int j, i, has_cells=0, has_lenses=0;

  /* finding the values of the current cell to decide if it is to be included in the calculation*/
  float cell_x = (float)(*tree)->center_mass_x/(float)(*tree)->total_mass;
  float cell_y = (float)(*tree)->center_mass_y/(float)(*tree)->total_mass;
  float dist = sqrt(pow(ray_x - cell_x,2) + pow(ray_y - cell_y, 2));
  float ratio = ((*tree)->width)/dist;
  //printf("ratio: %.8f\n", ratio);

  // checking if the current cell has any cells or lenses
  for(j=0; j<QUAD_IDX; j++){
    if((*tree)->lenses[j]) has_lenses++;
    if((*tree)->desc[j]) has_cells++;
  }

  if(delta == 0.0 && has_lenses>0){
    //lens_index += has_lenses;
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->lenses[i]){
        (*l_x)[lens_index] = (*tree)->lenses[i]->x;
        (*l_y)[lens_index] = (*tree)->lenses[i]->y;
        (*l_m)[lens_index] = (*tree)->lenses[i]->m;
        lens_index++;
      }
      if((*tree)->desc[i]) get_lenses(l_x, l_y, l_m, &((*tree)->desc[i]), delta, ray_x, ray_y);
    }
  }
  else if(has_lenses==1 && has_cells==0){
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->lenses[i]){
        (*l_x)[lens_index] = (*tree)->lenses[i]->x;
        (*l_y)[lens_index] = (*tree)->lenses[i]->y;
        (*l_m)[lens_index] = (*tree)->lenses[i]->m;
        lens_index++;
      }
    }
    return;
  }
  else if(has_lenses>0 && has_cells==0 && ratio>delta){
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->lenses[i]){
        (*l_x)[lens_index] = (*tree)->lenses[i]->x;
        (*l_y)[lens_index] = (*tree)->lenses[i]->y;
        (*l_m)[lens_index] = (*tree)->lenses[i]->m;
        lens_index++;
      }
    }
    return;
  }
  else if(has_lenses>0 && has_cells>0 && ratio>delta){
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->lenses[i]){
        (*l_x)[lens_index] = (*tree)->lenses[i]->x;
        (*l_y)[lens_index] = (*tree)->lenses[i]->y;
        (*l_m)[lens_index] = (*tree)->lenses[i]->m;
        lens_index++;
      }
      if((*tree)->desc[i]) get_lenses(l_x, l_y, l_m, &((*tree)->desc[i]), delta, ray_x, ray_y);
    }
    return;
  }
  else if(ratio<delta){
    (*l_x)[lens_index] = cell_x;
    (*l_y)[lens_index] = cell_y;
    (*l_m)[lens_index] =(*tree)->total_mass;
    lens_index++;
    return;
  }
  else{
    int i;
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->desc[i]) get_lenses(l_x, l_y, l_m, &((*tree)->desc[i]), delta, ray_x, ray_y);
    }
  }
}


__global__ void curand_setup(curandState* globalState, unsigned long seed){
	const unsigned int row = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int col = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int id = row*PIXEL_SIZE + col;
	curand_init(seed, id, 0, &globalState[id]);
}

__global__ void glensing(const float *d_l_x, const float *d_l_y, const float *d_l_mass, const size_t nobjects, unsigned int* results, 
							const d_constants* v, curandState* globalState, long int seed) {
	const unsigned int row = blockIdx.x*blockDim.x + threadIdx.x;
	const unsigned int col = blockIdx.y*blockDim.y + threadIdx.y;
	const unsigned int id = blockIdx.x*blockDim.x + blockDim.y;
	const float initial_x = (-v->image_scale_x) + row*v->increment_x;
	const float initial_y = (-v->image_scale_y) + col*v->increment_y;
	const float increment_x = v->increment_x;
	const float increment_y = v->increment_y;
	const float source_scale = v->source_scale;
	
	unsigned int lens_idx;
	if(IF_TREE_OPTIMIZATION)
		lens_idx = blockIdx.y*blockDim.y + blockIdx.x*nobjects; 
	else
		lens_idx = 0;
		
	float start_x, start_y, dx, dy;
	size_t k, it;
	__device__ __shared__ curandState localState;
	localState = globalState[id];
	
	for(it=0; it<RPP_Kernel; ++it){
		start_x = initial_x + curand_uniform(&localState) * increment_x;
		start_y = initial_y + curand_uniform(&localState) * increment_y;
		globalState[id] = localState;
		dx = (1-v->gamma_)*start_x - v->kappa_c*start_x;
		dy = (1+v->gamma_)*start_y - v->kappa_c*start_y;
	
		for(k = 0; k < nobjects; ++k) {
			float dist = pow(start_x - d_l_x[lens_idx + k], 2) + pow(start_y - d_l_y[lens_idx + k], 2);
			dx -= d_l_mass[lens_idx + k] * (start_x - d_l_x[lens_idx + k]) / dist;
			dy -= d_l_mass[lens_idx + k] * (start_y - d_l_y[lens_idx + k]) / dist;
		}
	
		if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
			(dy >= -source_scale/2) && (dy <= source_scale/2)) {
			int px = (dx + source_scale/2) / (source_scale/PIXEL_SIZE);
			int py = PIXEL_SIZE - (dy + source_scale/2) / (source_scale/PIXEL_SIZE);
			atomicAdd(&results[py * PIXEL_SIZE + px], 1);
		}
	}
}

int main(int argc, char** argv){	
	int i, j, m, array_size;
 	float increment_x, increment_y;
	root = NULL;
	lens *temp_lens;
	printf("PIXEL SIZE %d x %d\n", PIXEL_SIZE, PIXEL_SIZE);
	
// Load relevant settings and data
	if (argc < 2) error("Requires argument with lens positions and optional mass");
	setup_constants();
	
	d_constants *const_struct = (d_constants*)salloc(sizeof(d_constants));
	init_variables(const_struct);
	read_lenses(argv[1]);
	make_root(&root);

  /* creating the tree*/
  for(i=0; i < nobjects; i++){
    //temp = (cell *)salloc(sizeof(cell));                /*allocating memory to temp cell*/
    temp_lens = (lens *)salloc(sizeof(lens));                /*allocating memory to temp lens*/
    temp_lens->x = lens_x[i];                          /*assigning x-coord to cell_lens x-coord*/
    temp_lens->y = lens_y[i];                          /*assigning y-coord to cell_lens y-coord*/
    temp_lens->m = lens_mass[i];                          /*assigning mass*/
    make_tree(&root, temp_lens);
    free(temp_lens);
  }

  //print_tree(&root);

	printf("X %f and Y %f\n", image_scale_x, image_scale_y);
	increment_x = (image_scale_x * 2) / (PIXEL_SIZE);
	increment_y = (image_scale_y * 2) / (PIXEL_SIZE);
	const_struct->increment_x = increment_x;
	const_struct->increment_y = increment_y;
	printf("Increments for X %f and Y %f\n", increment_x, increment_y);
	
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
  printf("number of lenses: %d\n", *num_lenses);
  	lenses_x = (float *)salloc(sizeof(float)*array_size*(*num_lenses));
  	lenses_y = (float *)salloc(sizeof(float)*array_size*(*num_lenses));
  	lenses_m = (float *)salloc(sizeof(float)*array_size*(*num_lenses));
	
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
      		float *inc_lenses_x = (float *)salloc(sizeof(float)*(*num_lenses));
      		float *inc_lenses_y = (float *)salloc(sizeof(float)*(*num_lenses));
      		float *inc_lenses_m = (float *)salloc(sizeof(float)*(*num_lenses));
      		get_lenses(&inc_lenses_x, &inc_lenses_y, &inc_lenses_m, &root, delta, tile_x, tile_y);

      		for(m = 0; m<(*num_lenses); m++){
        		int idx = j*(PIXEL_SIZE/TILE_SIZE)*(*num_lenses) +i*(*num_lenses) +m;
        		lenses_x[idx] = inc_lenses_x[m];
        		lenses_y[idx] = inc_lenses_y[m];
        		lenses_m[idx] = inc_lenses_m[m];
      		}

      	free(inc_lenses_x);
      	free(inc_lenses_y);
      	free(inc_lenses_m);
      	count++;
    	}
  	}
	
	printf("size of lenses_x is %d\n", (int)sizeof(lenses_x));
	printf("nobjects is %d\n", (int)nobjects);
	printf("the size is %d\n", (*num_lenses)* GRID_SIZE * GRID_SIZE);
	
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
	
	if(nobjects != *num_lenses){
		printf("The number of lenses were reduced from %d to %d\n", (int)nobjects, *num_lenses);
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
		int error = 0;
				cudaMalloc(&d_const_struct, sizeof(d_constants));
		cudaMalloc(&globalState, sizeof(curandState)*GRID_SIZE*GRID_SIZE);
		
		dim3 bdim(1, 1);
		dim3 gdim(GRID_SIZE, GRID_SIZE); 
		curand_setup<<<gdim, bdim>>>(globalState, time(NULL));
		
		cudaMemcpy(d_const_struct, const_struct, sizeof(d_constants), cudaMemcpyHostToDevice);
			
		if(IF_TREE_OPTIMIZATION){
			nobjects = *num_lenses;
			cudaMalloc(&d_lens_x, sizeof(float) * nobjects * GRID_SIZE * GRID_SIZE);
			cudaMalloc(&d_lens_y, sizeof(float) * nobjects * GRID_SIZE * GRID_SIZE);
			cudaMalloc(&d_lens_mass, sizeof(float) * nobjects * GRID_SIZE * GRID_SIZE);
			cudaMemcpy(d_lens_x, lenses_x, sizeof(float) * nobjects * GRID_SIZE * GRID_SIZE, cudaMemcpyHostToDevice);
			cudaMemcpy(d_lens_y, lenses_y, sizeof(float) * nobjects * GRID_SIZE * GRID_SIZE, cudaMemcpyHostToDevice);
			cudaMemcpy(d_lens_mass, lenses_m, sizeof(float) * nobjects * GRID_SIZE * GRID_SIZE, cudaMemcpyHostToDevice);
		}
		else{
			cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
			cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
			cudaMalloc(&d_lens_mass, sizeof(float) * nobjects);
			cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
			cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
			cudaMemcpy(d_lens_mass, lens_mass, sizeof(float) * nobjects, cudaMemcpyHostToDevice);	
		}

		int iteration_No = 1;
		
		for(; iteration_No<=ITERATION&&!error; ++iteration_No){
			unsigned int *d_results;
			cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
			cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
			unsigned int *temp_results = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
			printf("OpenMP thread %d on Device %d; Iteration %d:%d \n", thread_No, device_No, iteration_No, ITERATION);
			dim3 bdim(TILE_SIZE, TILE_SIZE);
			dim3 gdim(GRID_SIZE, GRID_SIZE); 
			glensing<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_const_struct, globalState, time(NULL) + rand());
			cudaMemcpy(temp_results, d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost); 

			for(int k = 0; k<PIXEL_SIZE*PIXEL_SIZE; ++k){
				results[thread_No][k] += temp_results[k];
			}
			int total_x = total_r(temp_results, PIXEL_SIZE * PIXEL_SIZE);
			printf("OpenMP thread %d on Device %d; Iteration %d:%d Total Lenses observed %d\n", device_No, thread_No, iteration_No, ITERATION, total_x);
			free(temp_results);
			cudaFree(d_results);
		}
		
		
		/*free the variables*/
		cudaFree(d_lens_x);	  
		cudaFree(d_lens_y);
		cudaFree(d_lens_mass);
		cudaFree(d_const_struct);
		cudaFree(globalState);
	}	

	int r_c=0, t;
	for(r_c = 0; r_c < PIXEL_SIZE*PIXEL_SIZE; ++r_c){
		for(t=0; t<KERNEL_CALL_NUM; ++t)
			final_result[r_c] += results[t][r_c];
	}
	int total = total_r(final_result, PIXEL_SIZE * PIXEL_SIZE);
	printf("The total num of rays is %d \n", total);

	int highest_c = highest(final_result, PIXEL_SIZE * PIXEL_SIZE);
	write_pgm(final_result, PIXEL_SIZE, PIXEL_SIZE, highest_c);
  
	/*free variables*/
	free(lens_x);
	free(lens_y);
	free(lens_mass);
	free(const_struct);
	free(final_result);
	free(lenses_x);
	free(lenses_y);
	free(lenses_m);
	for(i=0; i<KERNEL_CALL_NUM; ++i)
		free(results[i]);

	cudaThreadExit();
	/*exit the program*/
	return 0;

}
