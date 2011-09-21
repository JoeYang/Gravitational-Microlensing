#include "utils.h"
#include "global.h"
#include "constants.h"

#define BLOCK_SIZE	(16)
#define GRID_SIZE	(16)
#define PIXEL_SIZE	(512)
#define show_lenses (1)

float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;

/*device lens_x and lens_y*/
float *d_lens_x;
float *d_lens_y;

void read_lenses(const char *filename) {
  size_t i, len = 0;
  char c, *line = NULL;
  FILE *fp;

  fprintf(stderr, "Reading in lenses...\n");
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

int highest(int *results, int size){
	int i = 0, highest_count=0;
	for(; i<size; ++i){
		if (results[i] > highest_count) 
			highest_count = results[i];
	}
	return highest_count;
}

__global__ void lens_position(const float *lens_x, const float *lens_y, int* results, vars* variables, int pixel_x, int pixel_y){
	const int x = threadIdx.x + blockIdx.x*blockDim.x;
  	const int y = threadIdx.y + blockIdx.y*blockDim.y;
  	const int i = y*BLOCK_SIZE*GRID_SIZE + x;
	const float dx = lens_x[i];
   	const float dy = lens_y[i];
   	float source_scale = variables->source_scale;
   	
  	if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
        (dy >= -source_scale/2) && (dy <= source_scale/2)) {
   	 	int px = (dx + source_scale/2) / (source_scale/pixel_x);
      	int py = (dy + source_scale/2) / (source_scale/pixel_y);
      	results[py * pixel_x + px] += 1;
  	}   
}

int main(int argc, char** argv)
{
	float increment_x, increment_y;
	int pixel_x = 512, pixel_y = 512;	
	// Load relevant settings and data
	if (argc < 2) error("Requires argument with lens positions and optional mass");
	vars *variables = setup_constants();
	read_lenses(argv[1]);	
	
	fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
	increment_x = (image_scale_x * 2) / pixel_x;
	increment_y = (image_scale_y * 2) / pixel_y;
	fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);
	
	int *results = (int *)calloc(pixel_x * pixel_y, sizeof(float));
	int *d_results;
	if (!results) error("calloc failed in allocating the result array");
	int highest_c = 0;
	
	// Fire off the light rays and record their end locations
 	cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
 	cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
 	cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));
 	
 	cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));	
  	cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
  	cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
      
  	dim3 dimb(BLOCK_SIZE, BLOCK_SIZE);
  	dim3 dimg_(GRID_SIZE, GRID_SIZE);
  	dim3 dimg(PIXEL_SIZE);
  	
  	//kernel	
  	// Place the gravitational objects in their relative positions
	if (show_lenses) {
		printf("oh yeah!\n");
	  	lens_position<<<dimg_, dimb>>>(d_lens_x, d_lens_y, d_results, variables, pixel_x, pixel_y);
	}
  	
  	cudaMemcpy(results, d_results,PIXEL_SIZE*PIXEL_SIZE*sizeof(float), cudaMemcpyDeviceToHost);

  	cudaFree(d_lens_x);
 	cudaFree(d_lens_y);
	cudaFree(d_results);
	
		
	highest_c = highest(results, pixel_x*pixel_y);
	assert(highest_c > 0 && "No pixels were written on the output map");
	write_pgm(results, pixel_x, pixel_y, highest_c);
	
	// Free the memory allocated during processing
	free(lens_x);
	free(lens_y);
	free(lens_mass);
	free(results);

  return 0;

}
