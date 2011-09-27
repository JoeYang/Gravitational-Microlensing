#include "global.h"

#ifndef UTIL_H_
#define UTIL_H_

void *salloc(size_t bytes);
void error(const char *msg);
void nextline(FILE *input, char **buf, size_t *len);

#endif
