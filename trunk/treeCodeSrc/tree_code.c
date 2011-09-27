#include "tree_code.h"

void setup_root(float, float, float, float, int);
void build_tree(int);
void calculate_cm(cell *, int);
void *salloc(size_t);
void print_tree(cell *, int);
void free_tree(cell *, int);
void remove_empty_cells(cell *, int);

cell * root;
cell * curr;
float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;
cell *cells; /* Cells to include in the lensing calculation */

int main(void) {
  int i;
  int ncells = 4;
  nobjects = 10;
  lens_x = (float *)salloc(sizeof(float) * nobjects);
  lens_y = (float *)salloc(sizeof(float) * nobjects);
  lens_mass = (float *)salloc(sizeof(float) * nobjects);
  
  for (i = 0; i < nobjects; ++i) {
    lens_x[i] = (i+1)*i;
    lens_y[i] = (i+1)*((float)i/2);
    lens_mass[i] = i+1;
  }
  
  setup_root(0, 1000, 1000, 0, ncells);
  build_tree(ncells);
  calculate_cm(root, ncells);
  
  printf("index, centre mass x, centre mass y, total mass\n"); 
  remove_empty_cells(root, ncells);
  print_tree(root, ncells);
  
  //free_tree(root);
  free(curr);
  
  return 0;
}

/* Remove cells that were created during the tree building process but were never assigned a lens */
void remove_empty_cells(cell *cellptr, int ncells) {
  int i;
  if (cellptr != 0) {
    for (i = 0; i < ncells; ++i) {
      if (cellptr->subcells[i] != 0) {
        if (cellptr->subcells[i]->mass == 0) {
          cellptr->subcells[i] = 0;
          free(cellptr->subcells[i]);
        }
        remove_empty_cells(cellptr->subcells[i], ncells);
      }
    }
  }
}

/* Print out the index, centre of mass and mass for each cell in the tree */
void print_tree(cell *cellptr, int ncells) {
  int i;
  for (i = 0; i < ncells; ++i) {
    if (cellptr->subcells[i] != 0) {
      print_tree(cellptr->subcells[i], ncells);
    }
  }
  printf("%d %f %f %f\n", cellptr->index, cellptr->cm_x, cellptr->cm_y, cellptr->mass);
}

/* Free all memory allocated to cells in the tree */
void free_tree(cell *cellptr, int ncells) {
  int i;
  if (cellptr != 0) {
    for (i = 0; i < ncells; ++i) {
      if (cellptr->subcells[i] != 0) {
        free_tree(cellptr->subcells[i], ncells);
      }
      if (cellptr->index != 0) {
        free(cellptr);
      }
    }
  }
}

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

/* Divide a cell by creating its 4 subcells */
void divide_cell(cell *curr, int ncells) {
	int i, j;
	for (j = 0; j < ncells; ++j) {
		cell * temp = (cell *)salloc(sizeof(cell));
		
		/* Calculate the index of the subcell */
		temp->index = ncells * curr->index + (j+1);
		
		/* Calculate the dimensions of the subcell */
		temp->top_left[0] = (curr->bottom_right[0] - curr->top_left[0])/2 * (j % 2) + curr->top_left[0];
		temp->top_left[1] = j > 1 ? (curr->top_left[1] - curr->bottom_right[1])/2 + curr->bottom_right[1] : curr->top_left[1];
		temp->bottom_right[0] = (j%2) == 0 ? (curr->bottom_right[0] - curr->top_left[0])/2 + curr->top_left[0] : curr->bottom_right[0];
		temp->bottom_right[1] = j < 2 ? (curr->top_left[1] - curr->bottom_right[1])/2 + curr->bottom_right[1] : curr->bottom_right[1];
		
		temp->cm_x = 0;
		temp->cm_y = 0;
		for (i = 0; i < 4; ++i) {
			temp->subcells[i] = 0;
		}
		
		curr->subcells[j] = temp;
	}
}

/* Finish calculating the centre of mass by dividing the weighted positions by the total mass */
void calculate_cm(cell *curr, int ncells) {
	int i;
	for (i = 0; i < ncells; ++i) {
		if (curr->subcells[i] > 0) {
			calculate_cm(curr->subcells[i], ncells);
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
		//include this cell/ lens in the calculation
	}
	// don't include this cell/ lens in the calculation
}

/* Construct the quadtree by adding all the lenses to it */
void build_tree(int ncells) {
	int i, j, lens1_index, lens2_index;
	cell *temp;
	
	for (i = 0; i < nobjects; ++i) {
		float dx = lens_x[i];
		float dy = lens_y[i];
		float lmass = lens_mass[i];
    
		/* Checks if lens is contained within the current cell */
		while (dx > curr->top_left[0] && dx < curr->bottom_right[0] &&
           dy > curr->bottom_right[1] && dy < curr->top_left[1]) {
			
			/* Update the total mass and weighted positions (to be used in centre of mass calculation) */
			curr->mass += lmass;
			curr->cm_y += lmass * dy;
			curr->cm_x += lmass * dx;
      
			/* If there are no lenses or subcells in the cell, place the lens in the cell. */
			if (!(curr->clens.pos_x) && !curr->subcells[0]) {
				
				curr->clens.pos_x = dx;
				curr->clens.pos_y = dy;
				curr->clens.mass = lmass;
				
				break;
        
			} else {
				/* More than one lens is being placed in a cell. If the cell has subcells, attempt to place the lens in one of those. Otherwise, divide the cell. */
				if (curr->subcells[0]) {
					
					for (j = 0; j < ncells; ++j) {
						temp = curr->subcells[j];
						if (dx > temp->top_left[0] && dx < temp->bottom_right[0] &&
								dy > temp->bottom_right[1] && dy < temp->top_left[1]) {
							
							curr = temp;
							break;
						}
					}
				} else {
					/* Divide a cell into subcells and assign lenses to their new subcells */
					do {
						lens1_index = 0;
						lens2_index = 0;
						
						divide_cell(curr, ncells);
						
						/* Put the two lenses in the cell into the subcells i.e. the original and new lens */
						for (j = 0; j < ncells; ++j) {
							temp = curr->subcells[j];
							
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
							curr = temp;
						}
						
					} while (lens1_index == lens2_index);
					
					break;
				}
			}
		}
		curr = root;
	}
}
