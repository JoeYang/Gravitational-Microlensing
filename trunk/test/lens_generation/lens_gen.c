#include "constants.h"
#include <math.h>

void lens_generate(float *new_lens_x, float *new_lens_y, int nobjects){
  
  setup_constants();

  int it = 0;
  srand(time(NULL));

  while(it < nobjects){
    float x, y;
    x = 2*lens_rad_x*rand() / RAND_MAX;
    y = 2*lens_rad_y*rand() / RAND_MAX;
    x -= lens_rad_x;
    y -= lens_rad_y;
    // Ensure the point is within the cirlce of radius lens_rad
    if (pow(x,2) + pow(y,2) <= pow(lens_rad,2)) {
      new_lens_x[it] = x;
      new_lens_y[it] = y;
      ++it;
    }
  }
}
