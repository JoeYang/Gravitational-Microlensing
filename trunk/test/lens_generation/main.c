#include <stdio.h>
#include "lens_gen.c"
#include "util.h"

float *lens_x, *lens_y, *lens_mass;
size_t nobjects;

int main(int argc, char *argv[]) {
  size_t it;

  if (argc < 2) {
    error("Please specify the number of lenses to generate");
  }
  nobjects = atoi(argv[1]);

  float *lens_x = (float*)salloc(nobjects*sizeof(float));
  float *lens_y = (float*)salloc(nobjects*sizeof(float));

  lens_generate(lens_x, lens_y, nobjects);

  for(it = 0; it < nobjects; it++){
    printf("%f %f\n", lens_x[it], lens_y[it]);
  }

  return 0;
}
