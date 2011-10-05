// All C includes must be wrapped in extern "C"
extern "C" {
#include "global.h"
#include "util.h"
#include "constants.h"
}
#include <assert.h>
#include <stdio.h>

#define PIXEL_SIZE (512)
#define TILE_SIZE (16)
#define GRID_SIZE (PIXEL_SIZE/TILE_SIZE)

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
  var->rpp = 32;
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

__global__ void group_glensing(const float *lens_x, const float *lens_y, const float *lens_mass, const size_t nobjects, unsigned int* results, const vars* v) {
  const unsigned int row = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int col = blockIdx.y*blockDim.y + threadIdx.y;
  
  const float initial_x = (-v->image_scale_x) + row*v->increment_x;
  const float initial_y = (-v->image_scale_y) + col*v->increment_y;

  const size_t ARRSIZE = 100;
  const unsigned int uniform_box = sqrtf((float)v->rpp);

  float start_x[ARRSIZE], start_y[ARRSIZE];
  float group_x[ARRSIZE], group_y[ARRSIZE];
  unsigned int it, noise_x, noise_y;
  size_t j, k;
  unsigned int px, py;
  float dist;

  for(it = 0; it < v->rpp/ARRSIZE; ++it) {
    // Set up the co-ordinates for the group
    for(j = 0; j < ARRSIZE; ++j) {
      noise_x = ((it * ARRSIZE) + j) % uniform_box;
      noise_y = ((it * ARRSIZE) + j) / uniform_box;
      start_x[j] = initial_x + noise_x * (v->increment_x / uniform_box);
      start_y[j] = initial_y + noise_y * (v->increment_y / uniform_box);
      group_x[j] = (1-v->gamma_)*start_x[j] - v->kappa_c*start_x[j];
      group_y[j] = (1+v->gamma_)*start_y[j] - v->kappa_c*start_y[j];
    }
    
    // Calculate the impact of each lens
    float lm, lx, ly;
    for(k = 0; k < nobjects; ++k) {
      lx = lens_x[k];
      ly = lens_y[k];
      lm = lens_mass[k];
      for(j = 0; j < ARRSIZE; ++j) {
        dist = pow(start_x[j] - lx, 2) + pow(start_y[j] - ly, 2);
        group_x[j] -= lm * (start_x[j] - lx) / dist;
        group_x[j] -= lm * (start_y[j] - ly) / dist;
      }
    }

    // Plot the output for each of the rays
    for(j = 0; j < ARRSIZE; ++j) {
      const float source_scale = v->source_scale;
      if ((group_x[j] >= -source_scale/2) && (group_x[j] <= source_scale/2) &&
          (group_y[j] >= -source_scale/2) && (group_y[j] <= source_scale/2)) {
        px = (group_x[j] + source_scale/2) / (source_scale/PIXEL_SIZE);
        py = PIXEL_SIZE - (group_y[j] + source_scale/2) / (source_scale/PIXEL_SIZE);
        atomicAdd(&results[py * PIXEL_SIZE + px], 1);
      }
    }

  }
}

__global__ void glensing(const float *lens_x, const float *lens_y, const float *lens_mass, const size_t nobjects, unsigned int* results, const vars* v) {
  const unsigned int row = blockIdx.x*blockDim.x + threadIdx.x;
  const unsigned int col = blockIdx.y*blockDim.y + threadIdx.y;
  
  const float initial_x = (-v->image_scale_x) + row*v->increment_x;
  const float initial_y = (-v->image_scale_y) + col*v->increment_y;

  const unsigned int uniform_box = sqrtf((float)v->rpp);

  float start_x, start_y, dx, dy;
  unsigned int it, noise_x, noise_y;
  size_t k;
  float dist;

  // TODO: Perform multiple ray calculations simultaneously
  // BUG: A larger value (> 100) of rpp results in a completely blank image
  for(it = 0; it < v->rpp; ++it) {
      noise_x = it % uniform_box;
      noise_y = it / uniform_box;
    start_x = initial_x + noise_x * v->increment_x / uniform_box;
    start_y = initial_y + noise_y * v->increment_y / uniform_box;

    dx = (1-v->gamma_)*start_x - v->kappa_c*start_x;
    dy = (1+v->gamma_)*start_y - v->kappa_c*start_y;

    for(k = 0; k < nobjects; ++k) {
      dist = pow(start_x - lens_x[k], 2) + pow(start_y - lens_y[k], 2);
      dx -= lens_mass[k] * (start_x - lens_x[k]) / dist;
      dy -= lens_mass[k] * (start_y - lens_y[k]) / dist;
    }

    const float source_scale = v->source_scale;
    if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
        (dy >= -source_scale/2) && (dy <= source_scale/2)) {
      int px = (dx + source_scale/2) / (source_scale/PIXEL_SIZE);
      int py = PIXEL_SIZE - (dy + source_scale/2) / (source_scale/PIXEL_SIZE);
      atomicAdd(&results[py * PIXEL_SIZE + px], 1);
      //results[py * PIXEL_SIZE + px] += 1;
    }
  }
}

int main(int argc, char** argv) {
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
  dim3 gdim(GRID_SIZE, GRID_SIZE);
  glensing<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_variables);
  //group_glensing<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_variables);
  cudaMemcpy(results, d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(unsigned int), cudaMemcpyDeviceToHost);

  int highest_c = highest(results, PIXEL_SIZE * PIXEL_SIZE);
  write_pgm(results, PIXEL_SIZE, PIXEL_SIZE, highest_c);

  // Free the memory allocated during processing
  // GPU
  cudaFree(d_lens_x);
  cudaFree(d_lens_y);
  cudaFree(d_lens_mass);
  cudaFree(d_results);
  cudaFree(d_variables);
  // CPU
  free(lens_x);
  free(lens_y);
  free(lens_mass);
  free(results);

  return 0;
}
