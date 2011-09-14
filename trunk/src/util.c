#include "util.h"

/* error • Reports an error and then exits with a failed status */
void error(const char *msg) {
  fprintf(stderr, "Error: %s\n", msg);
  exit(1);
}

/* salloc • Safe allocation of memory that exits cleanly in case of error */
void *salloc(size_t bytes) {
  void *ptr = malloc(bytes);
  /* Returning NULL for size 0 is valid */
  if (ptr == NULL && bytes != 0) {
    error("Cannot allocate object");
  }
  return ptr;
}

/* nextline • Loads next line of text, errors if no input */
void nextline(FILE *input, char **buf, size_t *len) {
  if (getline(buf, len, input) == -1) {
    error("getline has returned an error");
  }
}
