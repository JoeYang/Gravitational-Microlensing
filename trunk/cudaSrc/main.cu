// All C includes must be wrapped in extern "C"
extern "C" {
#include "global.h"
#include "util.h"
#include "constants.h"
}
#include <assert.h>
#include <stdio.h>

#define BLOCK_SIZE  (16)
#define SIZE    (16)
#define GRID_SIZE  (16)
#define PIXEL_SIZE  (512)
#define show_lenses (0)

float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;

/*device lens_x and lens_y*/
float *d_lens_x;
float *d_lens_y;

typedef struct vars {
  float kappa_star, kappa_c, gamma_, kappa, source_scale;
  float image_scale_fudge, image_scale_x, image_scale_y;
  float lens_scale_fudge_1, lens_scale_fudge_2, lens_rad_x, lens_rad_y, lens_rad;
  float lens_count;
  float increment_x, increment_y;
} vars;

void init_var(vars *var) {
  var->kappa_star = kappa_star;
  var->kappa_c = kappa_c;
  var->gamma_ = gamma_;
  var->kappa = kappa; 
  var->source_scale = source_scale;
  var->image_scale_fudge = image_scale_fudge;
  var->image_scale_x = image_scale_x; 
  var->image_scale_y = image_scale_y;
  var->lens_scale_fudge_1 = lens_scale_fudge_1; 
  var->lens_scale_fudge_2=lens_scale_fudge_2; 
  var->lens_rad_x = lens_rad_x;
  var->lens_rad_y = lens_rad_y;
  var->lens_rad = lens_rad;
  var->lens_count = lens_count;
  var->increment_x = 0;
  var->increment_y = 0;  
}

int highest(int *results, int size) {
  int i = 0, highest_count=0;
  for(; i<size; ++i){
    if (results[i] > highest_count) 
      highest_count = results[i];
  }
  return highest_count;
}

__global__ void glensing(const float *lens_x, const float *lens_y, int* results, int* row, vars* variables) {
  // x and y are source positions (in pixel)
  const int x = blockIdx.x; 
  const int y = *row;
  // i and j are light ray position on the lensing plane
  const int i = threadIdx.x;
  const int j = threadIdx.y;
  const float gamma_ = variables->gamma_;
  const float kappa_c = variables->kappa_c;
  const float image_scale_x = variables->image_scale_x;
  const float image_scale_y = variables->image_scale_y;
  const float source_scale = variables->source_scale;
  const float pixel_length_x = 2*image_scale_x/PIXEL_SIZE;
  const float pixel_length_y = 2*image_scale_y/PIXEL_SIZE;
  const float ray_dist_x = pixel_length_x/(2*BLOCK_SIZE);
  const float ray_dist_y = pixel_length_y/(2*BLOCK_SIZE);
  // dx and dy are on observation plane positions(in pixel)
  int dx = x, dy = y, iter;
  float dist;
  int b;
  for(b=0; b<4; b++){
    float start_x = (-1*image_scale_x) + x*pixel_length_x + i*ray_dist_x*(b%2+1);
    float start_y = (-1*image_scale_y) + y*pixel_length_y + j*ray_dist_y*(b/2+1);
     
    start_x *= (1-gamma_);
    start_y *= (1+gamma_);
    start_x -= kappa_c * start_x;
    start_y -= kappa_c * start_y;  
      
    for(iter = 0; iter < 2196; ++iter) {
      dist = pow(start_x - lens_x[iter], 2) + pow(start_y - lens_y[iter], 2);
      start_x -= (start_x - lens_x[iter]) / dist;
      start_y -= (start_y - lens_y[iter]) / dist;
    }
  
    if ((start_x >= -source_scale/2) && (start_x <= source_scale/2) &&
        (start_y >= -source_scale/2) && (start_y <= source_scale/2)) {     
      dx = (start_x + source_scale/2) / (source_scale/PIXEL_SIZE);
      dy = (start_y + source_scale/2) / (source_scale/PIXEL_SIZE);
      results[dy * PIXEL_SIZE + dx] += 1;
    }
  }  
}

int* init_row() {
  int *rows = (int*)malloc(PIXEL_SIZE*sizeof(int));
  int i;
  for(i=0; i<PIXEL_SIZE; i++) {
    rows[i] = i;
  }
  return rows;    
}

int main(int argc, char** argv) {
  float increment_x, increment_y;
  int pixel_x = 512, pixel_y = 512;  
  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  vars *variables = (vars *)salloc(sizeof(vars));
  init_var(variables);
  read_lenses(argv[1]);  
  
  fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
  increment_x = (image_scale_x * 2) / pixel_x;
  increment_y = (image_scale_y * 2) / pixel_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);
  
  int *results = (int *)calloc(pixel_x * pixel_y, sizeof(float));
  int *d_results;
  if (!results) error("calloc failed in allocating the result array");
  variables->increment_x = increment_x;
  variables->increment_y = increment_y;
  
  // cuda global memories
  vars *d_variables;
  int* row = init_row();
  int *d_row;
  cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
  cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
  cudaMalloc(&d_row, sizeof(int) * PIXEL_SIZE);
  cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));
  cudaMalloc(&d_variables, sizeof(vars));
   
  cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));  
  cudaMemcpy(d_variables, variables, sizeof(vars), cudaMemcpyHostToDevice);
  cudaMemcpy(d_row, row, PIXEL_SIZE*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
      
  dim3 dimb(SIZE, SIZE);
  dim3 dimt(GRID_SIZE, GRID_SIZE);
  dim3 dimg(PIXEL_SIZE);
    
  int r = 0;
  for(; r < PIXEL_SIZE; ++r) {
    glensing<<<dimg, dimb>>>(d_lens_x, d_lens_y, d_results, &d_row[r], d_variables);
  }
  cudaMemcpy(results, d_results,PIXEL_SIZE*PIXEL_SIZE*sizeof(float), cudaMemcpyDeviceToHost);
    
  int highest_c = highest(results, pixel_x * pixel_y);  
  write_pgm(results, pixel_x, pixel_y, highest_c);
  
  // Free the memory allocated during processing
  cudaFree(d_lens_x);
  cudaFree(d_lens_y);
  cudaFree(d_results);
  free(lens_x);
  free(lens_y);
  free(lens_mass);
  free(results);

  return 0;
}
