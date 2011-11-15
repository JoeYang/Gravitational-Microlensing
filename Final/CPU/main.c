#include "global.h"
#include "constants.h"
#include "util.h"

float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;

/* deflect • Deflect a given light ray according to the gravitational bodies along its path */
void deflect(float *x, float *y) {
  size_t i;
  float dist;
  // Start X / Y
  float start_x = *x, start_y = *y;
  *x = (1-gamma_)*start_x - kappa_c*start_x;
  *y = (1+gamma_)*start_y - kappa_c*start_y;
  for(i = 0; i < nobjects; ++i) {
    dist = pow(start_x - lens_x[i], 2) + pow(start_y - lens_y[i], 2);
    *x -= lens_mass[i] * (start_x - lens_x[i]) / dist;
    *y -= lens_mass[i] * (start_y - lens_y[i]) / dist;
  }
}

int main(int argc, const char *argv[])
{
  float increment_x, increment_y;
  int pixel_x = 512, pixel_y = 512, it;
  unsigned int rpp = 64;

  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  read_lenses(argv[1]);

  fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
  increment_x = (image_scale_x * 2) / pixel_x;
  increment_y = (image_scale_y * 2) / pixel_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);

  unsigned int *results = (int *)calloc(pixel_x * pixel_y, sizeof(unsigned int));
  if (!results) error("calloc failed in allocating the result array");
  int highest = 0;

  // Fire off the light rays and record their end locations
  int complete_iterations = 0;

  for(it = 0; it < rpp; ++it) {
    float x, y, dx, dy;
    const unsigned int uniform_box = sqrt(rpp);
    for(y = -image_scale_y; y < image_scale_y; y += increment_y) {
      for(x = -image_scale_x; x < image_scale_x; x += increment_x) {
        // Noise is uniformly distributed -- i.e. it's not really noise
        float noise_x = (it % uniform_box) * increment_x / uniform_box;
        float noise_y = (it / uniform_box) * increment_y / uniform_box;
        dx = x + noise_x;
        dy = y + noise_y;
        deflect(&dx, &dy);
        // Source plan (where collected) is -source_scale/2 to source_scale/2
        if ((dx > -source_scale/2) && (dx < source_scale/2) &&
            (dy > -source_scale/2) && (dy < source_scale/2)) {
          // Work out the nearest pixel to put this in to
          // Add to remove the negative part of source_scale and then source_scale / pixels
          int px = (dx + source_scale/2) / (source_scale/pixel_x);
          int py = pixel_y - (dy + source_scale/2) / (source_scale/pixel_y);
          results[py * pixel_x + px] += 1;
          if (results[py * pixel_x + px] > highest) highest = results[py * pixel_x + px];
        }
      }
    }
    fprintf(stderr, "\r%4d/%4d Complete\r", ++complete_iterations, rpp);
  }
  assert(highest > 0 && "No pixels were written on the output map");
  write_pgm(results, pixel_x, pixel_y, highest);

  // Free the memory allocated during processing
  free(lens_x);
  free(lens_y);
  free(lens_mass);
  free(results);

  return 0;
}
