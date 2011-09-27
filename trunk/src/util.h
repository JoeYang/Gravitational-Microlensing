#include "global.h"

#ifndef UTIL_H_
#define UTIL_H_

void *salloc(size_t bytes);
void error(const char *msg);
void nextline(FILE *input, char **buf, size_t *len);
void read_lenses(const char *filename);
void write_pgm(int *results, int pixel_x, int pixel_y, int highest);

#endif
