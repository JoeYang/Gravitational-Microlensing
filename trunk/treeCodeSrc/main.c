#include "global.h"
#include "constants.h"
#include "util.h"
#include "tree_code.h"

void initialise_array(int, int);
void calculate_cm(cell *, int);
void get_included_bodies(cell *, float, float, float);
void divide_cell(cell *, int);

cell ** cells;
float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;


/* read_lenses • Loads a lens file of format {x, y, (optional)mass} and allocate the correct sized array for the attributes */
void read_lenses(const char *filename) {
  size_t i, len = 0;
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
  float x, y, dx, dy, lmass, increment_x, increment_y;
  int pixel_x = 512, pixel_y = 512, it, i, j, index;
  cell * curr, * root, * temp;
  int ncells = 4; // number of subcells for a cell
  int lens1_index, lens2_index;

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
  //int show_lenses = 0;

  
	initialise_array(pixel_x, pixel_y);
	
	curr = cells[0];
	root = curr;
	
	for (i = 0; i < nobjects; ++i) {
		dx = lens_x[i];
		dy = lens_y[i];
		lmass = lens_mass[i];
		printf("Current lens: %f %f\n", dx, dy);
		
		/* Checks if lens is contained within the current cell */
		while (dx > curr->top_left[0] && dx < curr->bottom_right[0] &&
           dy > curr->bottom_right[1] && dy < curr->top_left[1]) {
			
			printf("Cell is within range\n");
			
			/* Update the total mass and weighted positions (to be used in centre of mass calculation) */
			curr->mass += lmass;
			curr->cm_x += lmass * dx;
			curr->cm_y += lmass * dy;
			
			/* If there are no lenses or subcells in the cell, place the lens in the cell. */
			if (!curr->clens.pos_x && !curr->subcells[0]) {
				printf("No lens in current cell\n");
				curr->clens.pos_x = dx;
				curr->clens.pos_y = dy;
				curr->clens.mass = lmass;
				
				break;
        
			} else {
				/* More than one lens is being placed in a cell. If the cell has subcells, attempt to place the lens in one of those. Otherwise, divide the cell. */
				if (curr->subcells[0]) {
					printf("Found subcells\n");
					
					for (j = 0; j < ncells; ++j) {
						temp = cells[curr->subcells[j]];
						if (dx > temp->top_left[0] && dx < temp->bottom_right[0] &&
							dy > temp->bottom_right[1] && dy < temp->top_left[1]) {
							
							curr = temp;
							break;
						}
					}
				} else {
					printf("No subcells. Need to create new cells\n");
					/* Divide a cell into subcells and assign lenses to their new subcells */
					do {
						lens1_index = 0;
						lens2_index = 0;
						
						divide_cell(curr, ncells);
						
						/* put the two lenses in the cell into the subcells i.e. the original and new lens */
						for (j = 0; j < ncells; ++j) {
							temp = cells[curr->subcells[j]];
							
							/* Find subcell for new lens */
							if (dx > temp->top_left[0] && dx < temp->bottom_right[0] &&
								dy > temp->bottom_right[1] && dy < temp->top_left[1]) {
								
								temp->clens.pos_x = dx;
								temp->clens.pos_y = dy;
								temp->clens.mass = lmass;
                
								temp->mass += lmass;
								temp->cm_x += lmass * dx;
								temp->cm_y += lmass * dy;
								
								lens1_index = temp->index;
							}
							
							/* Find subcell for lens originally in the cell */
							if (curr->clens.pos_x > temp->top_left[0] && curr->clens.pos_x < temp->bottom_right[0] &&
								curr->clens.pos_y > temp->bottom_right[1] && curr->clens.pos_y < temp->top_left[1]) {
                
								temp->clens = curr->clens;
								curr->clens.pos_x = 0;
								curr->clens.pos_y = 0;
								curr->clens.mass = 0;
								
								temp->mass += temp->clens.mass;
								temp->cm_x += temp->clens.mass * temp->clens.pos_x;
								temp->cm_y += temp->clens.mass * temp->clens.pos_y;
								
								lens2_index = temp->index;
							}
							
							if (lens1_index != 0 && lens2_index != 0) break;
						}
						
						if (lens1_index == lens2_index) {
							curr = cells[temp->index];
						}
						
					} while (lens1_index == lens2_index);
					
					break;
				}
			}
		}
		
		curr = root;
	}
	
	calculate_cm(root, ncells);
	/*
	printf("\nIndex of cell and its dimensions. i.e. top left corner and bottom right corner x, y co-ordinates\n");
	for (i = 0; i < 600; ++i) {
		if(cells[i]) {
			printf("%d %f %f %f %f\n", i, cells[i]->top_left[0], cells[i]->top_left[1], cells[i]->bottom_right[0], cells[i]->bottom_right[1]);
			printf("lens: %f %f %f\n", cells[i]->clens.pos_x, cells[i]->clens.pos_y, cells[i]->clens.mass);
			printf("cm: %f %f %f\n", cells[i]->mass, cells[i]->cm_x, cells[i]->cm_y);
		}
	}*/
  
  
  /*
  // Place the gravitational objects in their relative positions
  if (show_lenses) {
	  
	  
    for(it = 0; it < nobjects; ++it) {
      dx = lens_x[it];
      dy = lens_y[it];
      if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
          (dy >= -source_scale/2) && (dy <= source_scale/2)) {
        // Work out the nearest pixel to put this in to
        //result[py * pixel_x + px]
        // Add to remove the negative part of source_scale and then source_scale / pixels
        int px = (dx + source_scale/2) / (source_scale/pixel_x);
        int py = (dy + source_scale/2) / (source_scale/pixel_y);
        results[py * pixel_x + px] += 1;
        if (results[py * pixel_x + px] > highest) highest = results[py * pixel_x + px];
      }
    }
  }
*/
  // Fire off the light rays and record their end locations
  //#pragma omp parallel for
  for(it = 0; it < 50; ++it) {
    fprintf(stderr, "\nIteration Number %d Begun\n", it);
    for(y = -image_scale_y; y < image_scale_y; y += increment_y) {
      for(x = -image_scale_x; x < image_scale_x; x += increment_x) {
        // Noise should be up to one increment
        float noise = (((float)rand())/RAND_MAX * increment_x) - (increment_x/2);
        dx = x + noise;
        noise = (((float)rand())/RAND_MAX * increment_y) - (increment_y/2);
        dy = y + noise;
        //deflect(&dx, &dy);
        // Source plan (where collected) is -source_scale/2 to source_scale/2
        if ((dx >= -source_scale/2) && (dx <= source_scale/2) &&
            (dy >= -source_scale/2) && (dy <= source_scale/2)) {
          // Work out the nearest pixel to put this in to
          //result[py * pixel_x + px]
          // Add to remove the negative part of source_scale and then source_scale / pixels
          int px = (dx + source_scale/2) / (source_scale/pixel_x);
          int py = (dy + source_scale/2) / (source_scale/pixel_y);
          results[py * pixel_x + px] += 1;
          if (results[py * pixel_x + px] > highest) highest = results[py * pixel_x + px];
        }
      }
      // Provide an update counter in percent for our current progress
      // This doesn't play well with OpenMP for loops
      //fprintf(stderr, "\r%1.0f%% ", 100*(y+image_scale_y)/(image_scale_y*2));
    }
    fprintf(stderr, "Iteration Number %d Complete\n", it);
  }

  assert(highest > 0 && "No pixels were written on the output map");

  // Free the memory allocated during processing
  free(lens_x);
  free(lens_y);
  free(lens_mass);
  free(results);
  free(curr);

  return 0;
}

/* Create the array of cells and root node */
void initialise_array(int pixel_x, int pixel_y) {
	unsigned int max_index = 100; // TO DO: need to work out max number of indexes that can potentially be used
	cell * root;
	if ((cells = malloc(sizeof(cell) * max_index)) == NULL) {
		printf("Error allocating cells memory\n");
		exit(1);
	}

	if ((root = (cell *)malloc(sizeof(cell))) == NULL) {
		exit(1);
	}

	root->index = 0;
	root->top_left[0] = 0;
	root->top_left[1] = pixel_y;
	root->bottom_right[0] = pixel_x;
	root->bottom_right[1] = 0;
	root->mass = 0;
	root->cm_x = 0;
	root->cm_y = 0;
	
	cells[0] = root;
}

/* Divide a cell by creating its 4 subcells */
void divide_cell(cell *curr, int ncells) {
	int j;
	for (j = 0; j < ncells; ++j) {
		cell * temp = (cell *)malloc(sizeof(cell));
		
		/* calculate the index in which to store the subcell */
		temp->index = ncells * curr->index + (j+1);
		curr->subcells[j] = temp->index;
		
		/* calculate the dimensions of the subcell */
		temp->top_left[0] = (curr->bottom_right[0] - curr->top_left[0])/2 * (j % 2) + curr->top_left[0];
		temp->top_left[1] = j > 1 ? (curr->top_left[1] - curr->bottom_right[1])/2 + curr->bottom_right[1] : curr->top_left[1];
		temp->bottom_right[0] = (j%2) == 0 ? (curr->bottom_right[0] - curr->top_left[0])/2 + curr->top_left[0] : curr->bottom_right[0];
		temp->bottom_right[1] = j < 2 ? (curr->top_left[1] - curr->bottom_right[1])/2 + curr->bottom_right[1] : curr->bottom_right[1];
		
		//if (!(cells = realloc(cells, sizeof(cell) * sizeof(cells)+1))) {
		//	exit(1);
		//}
		
		cells[temp->index] = temp;
	}
}

/* Finish calculating the centre of mass by dividing the weighted positions by the total mass */
void calculate_cm(cell *curr, int ncells) {
	int i;
	for (i = 0; i < ncells; ++i) {
		if (curr->subcells[i] > 0) {
			calculate_cm(cells[curr->subcells[i]], ncells);
		}
	}
	if (curr->mass != 0) {
		curr->cm_x = curr->cm_x/curr->mass;
		curr->cm_y = curr->cm_y/curr->mass;
	}
}

/* Determine which bodies should be included in the deflection calculation */
void get_included_bodies(cell *curr, float accuracy, float ray_x, float ray_y) {
	float cell_width = curr->bottom_right[0] - curr->top_left[0];
	/* distance between the cell's centre of mass and the light ray */
	float d = sqrt(pow(ray_x - curr->cm_x, 2) + pow(ray_y - curr->cm_y, 2));
	
	/* if the ratio is less than the accuracy parameter, the body is sufficiently far away */
	if (cell_width/d > accuracy) {
		//include this cell/ lens calculation
	}
	// don't include this cell/ lens calculation
}
