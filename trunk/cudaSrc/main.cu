#include "utils.h"
#include "global.h"
#include "constants.h"
#include <assert.h>
#include <stdio.h>

#define BLOCK_SIZE	(16)
#define SIZE		(16)
#define GRID_SIZE	(16)
#define PIXEL_SIZE	(512)
#define show_lenses (0)

float *lens_x;
float *lens_y;
float *lens_mass;
int nobjects;

/*device lens_x and lens_y*/
float *d_lens_x;
float *d_lens_y;

void read_lenses(const char *filename) {
  size_t i;
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
  if (!(fout = fopen("img.pgm", "w"))) 
  	error("Can't open results file...");
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

__global__ void lens_position(const float *lens_x, const float *lens_y, int* results, vars* variables){
	const int column = threadIdx.x + blockIdx.x*blockDim.x;
  	const int row = threadIdx.y + blockIdx.y*blockDim.y;
  	const int i = row*SIZE*GRID_SIZE + column;
	const float source_scale = variables->source_scale;
	const int size = 2196;
	if(i<size){	
		float dx = lens_x[i];
	   	float dy = lens_y[i];
		if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
	 		(dy >= -source_scale/2) && (dy <= source_scale/2)) {  	  	
				int px = (dx + source_scale/2) / (source_scale/PIXEL_SIZE);
		 		int py = (dy + source_scale/2) / (source_scale/PIXEL_SIZE);
		 		results[py * PIXEL_SIZE + px] += 1;
		 	}
	}	   
}


__global__ void glensing(const float *lens_x, const float *lens_y, int* results, int* row, vars* variables){
  	// x and y are source positions (in pixel)
  	const int x = blockIdx.x; 
  	const int y = *row;
  	// i and j are light ray position on the lensing plane
  	const int i = threadIdx.x;
  	const int j = threadIdx.y;
  	const float gamma_ = variables->gamma_;
  	const float kappa_c = variables->kappa_c;
  	const float image_scale_x = variables->image_scale_x;
  	const float image_scale_y = variables->image_scale_y;
  	const float source_scale = variables->source_scale;
  	const float pixel_length_x = 2*image_scale_x/PIXEL_SIZE;
  	const float pixel_length_y = 2*image_scale_y/PIXEL_SIZE;
  	const float ray_dist_x = pixel_length_x/(2*BLOCK_SIZE);
  	const float ray_dist_y = pixel_length_y/(2*BLOCK_SIZE);
  	// dx and dy are on observation plane positions(in pixel)
  	int dx = x, dy = y, iter;
	float dist;
	int b;
 	for(b=0; b<4; b++){
	 	float start_x = (-1*image_scale_x) + x*pixel_length_x + i*ray_dist_x*(b%2+1);
	 	float start_y = (-1*image_scale_y) + y*pixel_length_y + j*ray_dist_y*(b/2+1);
	 	
	 	start_x *= (1-gamma_);
	  	start_y *= (1+gamma_);
	  	start_x -= kappa_c * start_x;
	  	start_y -= kappa_c * start_y;	
	  	
	  	for(iter = 0; iter < 2196; ++iter) {
	    	dist = pow(start_x - lens_x[iter], 2) + pow(start_y - lens_y[iter], 2);
		    start_x -= (start_x - lens_x[iter]) / dist;
	    	start_y -= (start_y - lens_y[iter]) / dist;
	 	}
	
	
		if ((start_x >= -source_scale/2) && (start_x <= source_scale/2) &&
	        (start_y >= -source_scale/2) && (start_y <= source_scale/2)) {     
	       dx = (start_x + source_scale/2) / (source_scale/PIXEL_SIZE);
	       dy = (start_y + source_scale/2) / (source_scale/PIXEL_SIZE);
	       results[dy * PIXEL_SIZE + dx] += 1;
			
		}
 	}	
}

int* init_row(){
	int *rows = (int*)malloc(PIXEL_SIZE*sizeof(int));
	int i;
	for(i=0; i<PIXEL_SIZE; i++)
		rows[i] = i;
	return rows;		
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
	variables->increment_x = increment_x;
	variables->increment_y = increment_y;
	
	// cuda global memories
	vars *d_variables;
	int* row = init_row();
	int *d_row;
 	cudaMalloc(&d_lens_x, sizeof(float) * nobjects);
 	cudaMalloc(&d_lens_y, sizeof(float) * nobjects);
 	cudaMalloc(&d_row, sizeof(int) * PIXEL_SIZE);
 	cudaMalloc(&d_results, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));
 	cudaMalloc(&d_variables, sizeof(vars));
 	
 	cudaMemset(d_results, 0, PIXEL_SIZE*PIXEL_SIZE*sizeof(float));	
  	cudaMemcpy(d_variables, variables, sizeof(vars), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_row, row, PIXEL_SIZE*sizeof(int), cudaMemcpyHostToDevice);
  	cudaMemcpy(d_lens_x, lens_x, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
  	cudaMemcpy(d_lens_y, lens_y, sizeof(float) * nobjects, cudaMemcpyHostToDevice);
      
  	dim3 dimb(SIZE, SIZE);
  	dim3 dimt(GRID_SIZE, GRID_SIZE);
  	dim3 dimg(PIXEL_SIZE);
  	
  	//kernel	
  	// Place the gravitational objects in their relative positions
	if (show_lenses) {
	  	lens_position<<<dimt, dimb>>>(d_lens_x, d_lens_y, d_results, d_variables);  	
	  	cudaThreadSynchronize();
	}
	else{
		int r = 0;
		for(; r < PIXEL_SIZE; ++r)
			glensing<<<dimg, dimb>>>(d_lens_x, d_lens_y, d_results, &d_row[r], d_variables);
	}  	
  	cudaMemcpy(results, d_results,PIXEL_SIZE*PIXEL_SIZE*sizeof(float), cudaMemcpyDeviceToHost);

  	cudaFree(d_lens_x);
 	cudaFree(d_lens_y);
	cudaFree(d_results);
	
		
	int highest_c = highest(results, pixel_x * pixel_y);
  	
	printf("the highest pixel count is %d\n", highest_c);
	
	write_pgm(results, pixel_x, pixel_y, highest_c);
	
	// Free the memory allocated during processing
	free(lens_x);
	free(lens_y);
	free(lens_mass);
	free(results);

  return 0;

}
