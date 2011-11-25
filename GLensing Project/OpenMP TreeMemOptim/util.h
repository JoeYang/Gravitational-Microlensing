#include "global.h"
#include "tree_struct.h"
#ifndef UTIL_H_
#define UTIL_H_

#include <time.h>
/*the struct that contains some key constant values that will be passed into GPU*/
typedef struct d_constants{
  unsigned int rpp;
  float kappa_c, gamma_, source_scale;
  float image_scale_x, image_scale_y;
  float increment_x, increment_y;
} d_constants;

void *salloc(size_t bytes);
void error(const char *msg);
void nextline(FILE *input, char **buf, size_t *len);
void read_lenses(const char *filename, cell ** root);
int highest(unsigned int *results, unsigned int size);
int total_r(unsigned int *results, unsigned int size);
void write_pgm(unsigned int *results, int pixel_x, int pixel_y, int highest);

#endif
