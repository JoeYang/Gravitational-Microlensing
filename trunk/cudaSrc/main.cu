#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>

#include "utils.h"
#include "constants.h"

#define BLOCK_SIZE	16
#define PIXEL_SIZE	512


#define  

float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;

/*device lens_x and lens_y*/
float *d_lens_x;
float *d_lens_y;

/* read_lenses • Loads a lens file of format {x, y, (optional)mass} and allocate the correct sized array for the attributes */
void read_lenses(const char *filename) {
  size_t i, len = 0;
  float x, y;
  char c, *tmp, *line = NULL;
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
    if(fscanf(fp, "%f %f", &x, &y)!=2) error("wrong data!\n");
    lens_x[i] = x;
    lens_y[i] = y;
    lens_mass[i] = 1;
    
  }

  if (fclose(fp) != 0) error("Can't close lens file...");
  // Deallocate memory used by line
  free(line);
} 

int highest(int *result, int size){
	int i = 0, highest_count=0;
	for(; i<size; ++i){
		if (results[i] > highest_count) 
			highest_count = results[i];
	}
	return highest_count;
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


__global__ void kernel(const float *lens_x, const float *lens_y, int* result, const int iter){
  	const int x = blockIdx.x; 
  	const int y = iter;
  	const int i = threadIdx.x;
  	const int j = threadIdx.y;
  	int dx = x, dy = y;
  	
  	/* deflection calculation*/
  
  	
  	/*finish calculation*/
  	
  	result[dy*PIXEL_SIZE + dx] ++;
}


int main(int argc, char** argv)
{
	float x, y, dx, dy, increment_x, increment_y;
	int pixel_x = 512, pixel_y = 512, it;
	
	// Load relevant settings and data
	if (argc < 2) error("Requires argument with lens positions and optional mass");
	setup_constants();
	read_lenses(argv[1]);
	fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
	increment_x = (image_scale_x * 2) / (pixel_x*10);
	increment_y = (image_scale_y * 2) / (pixel_y*10);
	fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);
	
	int *results = (int *)calloc(pixel_x * pixel_y, sizeof(float));
	int *d_results;
	if (!results) error("calloc failed in allocating the result array");
  
    /* here starts the cuda*/
    cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
 	cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
 	cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));
 	
 	cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));
  	
  	cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
  	cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
      
  	dim3 dimb(BLOCK_SIZE, BLOCK_SIZE);
  	dim3 dimg(PIXEL_SIZE);
  	
  	for(int i=0; i<=PIXEL_SIZE; ++i){
  		kernel<<<dimg, dimb>>>(d_lens_x, d_lens_y, d_results, i);
  		fprintf(stderr, "\r%1.0f%% ", 100*i*1.0/PIXEL_SIZE);
  	}	
  	
  	cudaMemcpy(d_results, results,PIXEL_SIZE*PIXEL_SIZE*sizeof(float), cudaMemcpyDeviceToHost);

  	cudaFree(d_lens_x);
 	cudaFree(d_lens_y);
	cudaFree(d_results);

  	
  	/* cuda kernel finishes*/
  	fprintf(stderr, "\n");
 	write_pgm(results, pixel_x, pixel_y, highest(results, PIXEL_SIZE*PIXEL_SIZE));  	
 	
	return 0;
}
