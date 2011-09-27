#include "global.h"
#include "constants.h"
#include "util.h"

float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;

/* read_lenses • Loads a lens file of format {x, y, (optional)mass} and allocate the correct sized array for the attributes */
void read_lenses(const char *filename) {
  size_t i, len = 0;
  char c, *tmp, *line = NULL;
  FILE *fp;

  fprintf(stderr, "Reading in xses...\n");
  if (!(fp = fopen(filename, "r"))) error("Can't open lens file...");

  nobjects = 0;
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
 	if(fscanf(fp, "%f %f", &lens_x[i], &lens_y[i])!=2)
    	error("invalid input!");
    lens_mass[i] = 1;	
  }

  if (fclose(fp) != 0) error("Can't close lens file...");
  // Deallocate memory used by line
  free(line);
}

/* get_pixel • Return the correct pixel in the output map */
// TODO

/* deflect • Deflect a given light ray according to the gravitational bodies along its path */
void deflect(float *x, float *y) {
  size_t i;
  float dist;
  // Start X / Y
  float start_x = *x, start_y = *y;
  *x *= (1-gamma_);
  *y *= (1+gamma_);
  *x -= kappa_c * start_x;
  *y -= kappa_c * start_y;
  for(i = 0; i < nobjects; ++i) {
    dist = pow(start_x - lens_x[i], 2) + pow(start_y - lens_y[i], 2);
    *x -= lens_mass[i] * (start_x - lens_x[i]) / dist;
    *y -= lens_mass[i] * (start_y - lens_y[i]) / dist;
  }
}

/* write_pgm • Output the results as a PGM (portable gray map) image for review */
void write_pgm(int *results, int pixel_x, int pixel_y, int highest) {
  FILE *fout;
  fprintf(stderr, "Writing resulting image...\n");
  if (!(fout = fopen("img.pgm", "w"))) error("Can't open results file...");
  // Writing the PGM format which starts with P2
  fprintf(fout, "P2\n");
  // Followed by pixel width, height and the value considered white
  fprintf(fout, "%d %d\n", pixel_x, pixel_y);
  fprintf(fout, "%d\n", highest);
  // Print each value in a row of WIDTH length
  int px, py;
  for(py = 0; py < pixel_y; ++py) {
    for(px = 0; px < pixel_x; ++px) {
      fprintf(fout, "%d ", results[py * pixel_x + px]);
    }
    fprintf(fout, "\n");
  }
  if (fclose(fout) != 0) error("Can't close results file...");
}

int main(int argc, const char *argv[])
{
  float increment_x, increment_y;
  int pixel_x = 512, pixel_y = 512, it;
  unsigned int rpp = 500;

  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  read_lenses(argv[1]);

  fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
  increment_x = (image_scale_x * 2) / pixel_x;
  increment_y = (image_scale_y * 2) / pixel_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);

  int *results = (int *)calloc(pixel_x * pixel_y, sizeof(float));
  if (!results) error("calloc failed in allocating the result array");
  int highest = 0;

  // Fire off the light rays and record their end locations
  int complete_iterations = 0;
  #pragma omp parallel for
  for(it = 0; it < rpp; ++it) {
    float x, y, dx, dy;
    for(y = -image_scale_y; y < image_scale_y; y += increment_y) {
      for(x = -image_scale_x; x < image_scale_x; x += increment_x) {
        // Noise is uniformly distributed -- i.e. it's not really noise
        float noise_x = it * increment_x / rpp;
        float noise_y = it * increment_y / rpp;
        dx = x + noise_x;
        dy = y + noise_y;
        deflect(&dx, &dy);
        // Source plan (where collected) is -source_scale/2 to source_scale/2
        if ((dx > -source_scale/2) && (dx < source_scale/2) &&
            (dy > -source_scale/2) && (dy < source_scale/2)) {
          // Work out the nearest pixel to put this in to
          // Add to remove the negative part of source_scale and then source_scale / pixels
          int px = (dx + source_scale/2) / (source_scale/pixel_x);
          int py = source_scale/ (source_scale/pixel_y) - (dy + source_scale/2) / (source_scale/pixel_y);
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
