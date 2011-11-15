// All C includes must be wrapped in extern "C"
extern "C" {
#include "global.h"
#include "util.h"
#include "constants.h"
#include "tree_struct.h"
}
#include "treecode_kernel.h"

#define pixel_x (512)
#define pixel_y (512)
#define QUAD_IDX (4)
#define PIXEL_SIZE (512)
#define TILE_SIZE (16)
#define GRID_SIZE (PIXEL_SIZE/TILE_SIZE)

float4 *lenses, *d_lenses;

float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;
cell *root; 	/* global pointer to root cell so we can get to it from anywhere */
const float delta = 0.0;		/*variable to determine accuracy*/
static int lens_index = 0;

void init_var(vars *var) {
  var->rpp = 50;
  var->kappa_c = kappa_c;
  var->gamma_ = gamma_;
  var->source_scale = source_scale;
  var->image_scale_x = image_scale_x;
  var->image_scale_y = image_scale_y;
  var->increment_x = 0;
  var->increment_y = 0;
}

/* Method to get the highest pixel count int he result set to pass to
* the pgm creater */
/*int highest(unsigned int *results, unsigned int size) {
  unsigned int i, highest_count = 0;
  for(i = 0; i < size; ++i){
    if (results[i] > highest_count)
      highest_count = results[i];
  }
  return highest_count;
}*/

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

/*
* main method
*/
int main(int argc, const char *argv[]){
  int i, j, m, array_size;
  float increment_x, increment_y, elapsedTime;
  root = NULL;
  lens *temp_lens;

  //variables to capture the start and stop time
  cudaEvent_t start, stop;

  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  vars *variables = (vars *)salloc(sizeof(vars));
  init_var(variables);
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
  
  fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
  increment_x = (image_scale_x * 2) / pixel_x;
  increment_y = (image_scale_y * 2) / pixel_y;
  variables->increment_x = increment_x;
  variables->increment_y = increment_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);

  /* setting up the array of lenses denoted by the accuracy chosen*/

  unsigned int *results = (unsigned int *)calloc(PIXEL_SIZE * PIXEL_SIZE, sizeof(unsigned int));
  //unsigned int *temp = (unsigned int *)salloc(sizeof(unsigned int)*TILE_SIZE*TILE_SIZE);
  unsigned int *d_results;
  if (!results) error("calloc failed in allocating the result array");

  int *num_lenses = (int *)salloc(sizeof(int));
  if((PIXEL_SIZE % TILE_SIZE) == 0){
    array_size = (PIXEL_SIZE/TILE_SIZE)*(PIXEL_SIZE/TILE_SIZE);
  }
  else{
    error("PIXEL_SIZE is not a multiple of TILE_SIZE");
  }

  //here we get an initial value of x and y so we can get the number of lenses here to allocate the array
  // that will be passed to the GPU
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

//moving all the cuda code here as we need to calculate all the walks before passing them all off to the GPU
// we can do this as the total number of interacting lenses with rays will still be far less than the
// total number of lenses read in

  //starting timer
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  // Setting up CUDA global memory
  vars *d_variables;
  cudaMalloc(&d_variables, sizeof(vars));
  cudaMemcpy(d_variables, variables, sizeof(vars), cudaMemcpyHostToDevice);
  cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));
  cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int));

  /* allocating space on the gpu for the array of lenses*/
  cudaMalloc(&d_lenses, sizeof(float4)*(*num_lenses)*(PIXEL_SIZE/TILE_SIZE)*(PIXEL_SIZE/TILE_SIZE));
  cudaMemset(d_lenses, 0, sizeof(float4)*(*num_lenses)*(PIXEL_SIZE/TILE_SIZE)*(PIXEL_SIZE/TILE_SIZE));
  cudaMemcpy(d_lenses, lenses, sizeof(float4)*(*num_lenses)*(PIXEL_SIZE/TILE_SIZE)*(PIXEL_SIZE/TILE_SIZE), cudaMemcpyHostToDevice);

  // Perform gravitational microlensing
  dim3 bdim(TILE_SIZE, TILE_SIZE);
  dim3 gdim(GRID_SIZE, GRID_SIZE);
  glensing<<<gdim, bdim>>>(d_lenses, (*num_lenses), d_results, d_variables);
  //group_glensing<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_variables);
  cudaMemcpy(results, d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost);

  //stopping timer
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);

  //getting elapsed time
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf("Total time taken on the GPU is: %f ms\n", elapsedTime);

  //destroying timers
  cudaEventDestroy(start);
  cudaEventDestroy(stop);

  int highest_c = highest(results, PIXEL_SIZE * PIXEL_SIZE);
  write_pgm(results, PIXEL_SIZE, PIXEL_SIZE, highest_c);

  // Free the memory allocated during processing
  //GPU
  cudaFree(d_variables);
  cudaFree(d_lenses);
  cudaFree(d_results);
  // CPU
  free(lenses);
  free(num_lenses);
  free(variables);
  free(lens_x);
  free(lens_y);
  free(lens_mass);
  //free(temp);
  free(results);
  free_tree(&root);
  return 0;
}
