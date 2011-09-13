#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "util.h"

/* read_lenses • Loads a lens file of format {x, y, (optional)mass} and allocate the correct sized array for the attributes */
void read_lenses(float *lens_x, float *lens_y, float *lens_mass, const char *filename) {
  size_t i, len = 0, nobjects = 0;
  char c, *tmp, *line = NULL;
  FILE *fp;

  fprintf(stderr, "Reading in lenses...\n");
  if (!(fp = fopen(filename, "r"))) error("Can't open lens file...");

  // Count the number of lenses we must allocate for (one per line)
  while ((c = getc(fp)) != EOF) {
    if (c == '\n') ++nobjects;
  }
  fprintf(stderr, "Total lenses found: %d\n", nobjects);
  // Seek to the start of the file for actual reading
  fseek(fp, 0, SEEK_SET);

  // Allocate memory for the lenses
  lens_x = (float *)salloc(sizeof(float) * nobjects);
  lens_y = (float *)salloc(sizeof(float) * nobjects);
  lens_mass = (float *)salloc(sizeof(float) * nobjects);

  for(i = 0; i < nobjects; ++i) {
    // Read object definition
    nextline(fp, &line, &len);
    // Object x and y values are required
    if ((tmp = (char *)strtok(line, " \n")) == NULL) error("Lens has no x");
    lens_x[i] = atof(tmp);
    if ((tmp = (char *)strtok(NULL, " \n")) == NULL) error("Lens has no y");
    lens_y[i] = atof(tmp);
    // Mass is optional
    if ((tmp = (char *)strtok(NULL, " \n")) == NULL) lens_mass[i] = 1;
    else lens_mass[i] = atof(tmp);
  }

  if (fclose(fp) != 0) error("Can't close lens file...");
  // Deallocate memory used by line
  free(line);
}

int main(int argc, const char *argv[])
{
  float *lens_x;
  float *lens_y;
  float *lens_mass;

  if (argc < 2) error("Requires argument with lens positions and optional mass");

  setup_constants();
  read_lenses(lens_x, lens_y, lens_mass, argv[1]);

  return 0;
}
