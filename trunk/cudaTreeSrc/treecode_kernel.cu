#include "global.h"
#include "util.h"
#include "constants.h"
#include "tree_struct.h"
#include "treecode_kernel.h"
#define PIXEL_SIZE (512)

__global__ void glensing(const float4 * lenses, const size_t nobjects, unsigned int* results, const vars* v) {
  const unsigned int lens_idx = blockIdx.y*blockDim.y + blockIdx.x*nobjects;  //+ blockIdx.y*blockDim.y;

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
      dist = pow(start_x - lenses[lens_idx + k].x, 2) + pow(start_y - lenses[lens_idx + k].y, 2);
      dx -= lenses[lens_idx + k].w * (start_x - lenses[lens_idx + k].x) / dist;
      dy -= lenses[lens_idx + k].w * (start_y - lenses[lens_idx + k].y) / dist;
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

