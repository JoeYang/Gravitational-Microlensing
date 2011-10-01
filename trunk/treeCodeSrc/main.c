#include "global.h"
#include "constants.h"
#include "util.h"
#include "tree_code.h"

float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;
cell *root;
cell *curr;
cell *cells;
int count;

/* deflect • Deflect a given light ray according to the gravitational bodies along its path */
void deflect(float *x, float *y) {
  size_t i;
  float dist;
  // Start X / Y
  float start_x = *x, start_y = *y;
  *x = (1-gamma_)*start_x - kappa_c*start_x;
  *y = (1+gamma_)*start_y - kappa_c*start_y;
  for(i = 0; i < count; ++i) {
    dist = pow(start_x - cells[i].cm_x, 2) + pow(start_y - cells[i].cm_y, 2);
    *x -= cells[i].mass * (start_x - cells[i].cm_x) / dist;
    *y -= cells[i].mass * (start_y - cells[i].cm_y) / dist;
  }
}

/* Set up the root node */ 
void setup_root(float top_x, float top_y, float bottom_x, float bottom_y, int ncells) {
  int i;
  curr = (cell *)salloc(sizeof(cell));
  curr->index = 0;
  curr->top_left[0] = top_x;
  curr->top_left[1] = top_y;
  curr->bottom_right[0] = bottom_x;
  curr->bottom_right[1] = bottom_y;
  curr->mass = 0;
  curr->cm_x = 0;
  curr->cm_y = 0;
  
  for (i = 0; i < ncells; ++i) {
    curr->subcells[i] = 0;
  }
  
  root = curr;
}

/* Count the number of bodies that should be included in the deflection calculation */
void count_included_bodies(cell *cellptr, float accuracy, float ray_x, float ray_y) {
  int i;
  int got_subcells = 0;
  for (i = 0; i < 4; ++i) {
  	if (cellptr->subcells[i] != 0) {
 	  got_subcells = 1;
	}
  }

  if (got_subcells) {
  	float cell_width = cellptr->bottom_right[0] - cellptr->top_left[0];
  	float d = sqrt(pow(ray_x - cellptr->cm_x, 2) + pow(ray_y - cellptr->cm_y, 2));
  	if ((cell_width/d) < accuracy) {
      count++;
    }
    else {
  	  for (i = 0; i < 4; ++i) {
  	    if (cellptr->subcells[i] != 0) {
  	      count_included_bodies(cellptr->subcells[i], accuracy, ray_x, ray_y);
	    }
	  }
    }
  }
  else {
    count++;
  }
}


/* Determine which bodies should be included in the deflection calculation */
void get_included_bodies(cell *cellptr, float accuracy, float ray_x, float ray_y) {
  int i;
  int got_subcells = 0;
  for (i = 0; i < 4; ++i) {
  	if (cellptr->subcells[i] != 0) {
 	  got_subcells = 1;
	}
  }

  if (got_subcells) {
  	float cell_width = cellptr->bottom_right[0] - cellptr->top_left[0];
  	float d = sqrt(pow(ray_x - cellptr->cm_x, 2) + pow(ray_y - cellptr->cm_y, 2));
  	if ((cell_width/d) < accuracy) {
      //increase_array(&cells, count);
      cells[count++] = *cellptr;
    }
    else {
  	  for (i = 0; i < 4; ++i) {
  	    if (cellptr->subcells[i] != 0) {
  	      get_included_bodies(cellptr->subcells[i], accuracy, ray_x, ray_y);
	    }
	  }
    }
  }
  else {
    //increase_array(&cells, count);
    cells[count++] = *cellptr;
  }
}

int main(int argc, const char *argv[])
{
  float increment_x, increment_y;
  int pixel_x = 512, pixel_y = 512, it;
  unsigned int rpp = 500;
  int ncells = 4;
  float accuracy = 0.6;
  count = 0;
  root = NULL;
  curr = NULL;
  
  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  read_lenses(argv[1]);
  
  setup_root(-lens_rad, lens_rad, lens_rad, -lens_rad, ncells);
  
  build_tree(ncells, root);
  calculate_cm(root, ncells);

  //printf("All cells\n");
  //printf("index, centre mass x, centre mass y, total mass\n"); 
  remove_empty_cells(root, ncells);
  //print_tree(root, ncells);

  fprintf(stderr, "X %f and Y %f\n", image_scale_x, image_scale_y);
  increment_x = (image_scale_x * 2) / pixel_x;
  increment_y = (image_scale_y * 2) / pixel_y;
  fprintf(stderr, "Increments for X %f and Y %f\n", increment_x, increment_y);

  unsigned int *results = (int *)calloc(pixel_x * pixel_y, sizeof(unsigned int));
  if (!results) error("calloc failed in allocating the result array");
  int highest = 0;

  // Fire off the light rays and record their end locations
  int complete_iterations = 0;
  //#pragma omp parallel for
  for(it = 0; it < rpp; ++it) {
    float x, y, dx, dy;
    for(y = -image_scale_y; y < image_scale_y; y += increment_y) {
      for(x = -image_scale_x; x < image_scale_x; x += increment_x) {
        // Noise is uniformly distributed -- i.e. it's not really noise
        float noise_x = it * increment_x / rpp;
        float noise_y = it * increment_y / rpp;
        dx = x + noise_x;
        dy = y + noise_y;
        
        count_included_bodies(root, accuracy, dx, dy);
        cells = (cell *)malloc(sizeof(cell) * count);
        count = 0;
        get_included_bodies(root, accuracy, dx, dy);
        
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
        count = 0;
        free(cells);
		cells = NULL;
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
  free(curr);

  return 0;
}
