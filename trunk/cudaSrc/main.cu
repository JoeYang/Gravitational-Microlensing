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
  var->rpp = 500;
  var->kappa_c = kappa_c;
  var->gamma_ = gamma_;
  var->source_scale = source_scale;
  var->image_scale_x = image_scale_x;
  var->image_scale_y = image_scale_y;
  var->increment_x = 0;
  var->increment_y = 0;
}

int highest(int *results, int size) {
  int i, highest_count = 0;
  for(i = 0; i < size; ++i){
    if (results[i] > highest_count)
      highest_count = results[i];
  }
  return highest_count;
}

__global__ void glensing(const float *lens_x, const float *lens_y, const float *lens_mass, const size_t nobjects, int* results, vars* v) {
  // Position of Block -- likely useful in tree calculations
  //const int bx = blockIdx.x;
  //const int by = blockIdx.y;
  //
  const unsigned int row = blockIdx.y*blockDim.y + threadIdx.y;
  const unsigned int col = blockIdx.x*blockDim.x + threadIdx.x;

  float start_x, start_y;
  float noise_x = v->increment_x / v->rpp;
  float noise_y = v->increment_y / v->rpp;
  int dx, dy, it;
  float dist;

  // TODO: Perform multiple ray calculations simultaneously
  for(it = 0; it < v->rpp; ++it) {
    start_x = (-v->image_scale_x) + row*v->increment_x;
    start_x += it * noise_x;
    start_y = (-v->image_scale_y) + col*v->increment_y;
    start_y += it * noise_y;
    dx = start_x;
    dy = start_y;

    dx = (1-v->gamma_)*start_x - v->kappa_c*start_x;
    dx = (1+v->gamma_)*start_y - v->kappa_c*start_y;

    for(iter = 0; iter < nobjects; ++iter) {
      dist = pow(start_x - lens_x[iter], 2) + pow(start_y - lens_y[iter], 2);
      start_x -= lens_mass[iter] * (start_x - lens_x[iter]) / dist;
      start_y -= lens_mass[iter] * (start_y - lens_y[iter]) / dist;
    }

    const float source_scale = v->source_scale;
    if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
        (dy >= -source_scale/2) && (dy <= source_scale/2)) {
      results[dy * PIXEL_SIZE + dx] += 1;
    }
  }
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
  variables->increment_x = increment_x;
  variables->increment_y = increment_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);

  int *results = (int *)calloc(pixel_x * pixel_y, sizeof(float));
  int *d_results;
  if (!results) error("calloc failed in allocating the result array");

  // Setting up CUDA global memory
  vars *d_variables;
  cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
  cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
  cudaMalloc(&d_lens_mass, sizeof(float) * nobjects);

  cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));
  cudaMalloc(&d_variables, sizeof(vars));

  cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));
  cudaMemcpy(d_variables, variables, sizeof(vars), cudaMemcpyHostToDevice);
  cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
  cudaMemcpy(d_lens_mass, lens_mass, sizeof(float) * nobjects, cudaMemcpyHostToDevice);

  // Perform gravitational microlensing
  dim3 bdim(TILE_SIZE, TILE_SIZE);
  dim3 gdim(GRID_SIZE, GRID_SIZE);
  glensing<<<gdim, bdim>>>(d_lens_x, d_lens_y, d_lens_mass, nobjects, d_results, d_variables);
  cudaMemcpy(results, d_results,PIXEL_SIZE*PIXEL_SIZE*sizeof(float), cudaMemcpyDeviceToHost);

  int highest_c = highest(results, pixel_x * pixel_y);
  write_pgm(results, pixel_x, pixel_y, highest_c);

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
