#include "tree_code.h"
#include "util.h"

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


void increase_array(cell **cell_array, int num_elements) {
  cell *temp = NULL;
  if ((temp = (cell *)realloc(*cell_array, (num_elements+1) * sizeof(cell))) == NULL) {
  	  printf("ERROR: realloc failed");
  }

  *cell_array = temp;
}

/* Construct the quadtree by adding all the lenses to it */
void build_tree(int ncells, cell * root) {
  
  float dx, dy, lmass;
  int i, j, lens1_index, lens2_index;
  cell *temp = NULL, *curr = root;

  for (i = 0; i < nobjects; ++i) {
	dx = lens_x[i];
	dy = lens_y[i];
	lmass = lens_mass[i];
   
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
