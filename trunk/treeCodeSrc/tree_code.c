#include "tree_code.h"

void calculate_r(cell **);

int main(void) {
	cell * curr, * root, * temp;
	int i, j, index;
	float dx, dy;
	float lmass;
	int lens1_index, lens2_index;
	int max_index = 200; // TO DO: need to work out max number of indexes that can potentially be used
	int ncells = 4; // number of cells to split a cell into when two lenses are in the same cell
	int nobjects = 5;
	
	/* test data */
	int pixel_x = 5, pixel_y = 5; // dimensions of the lensing plane
	float lens_x[5] = {4, 3.5, 4.5, 2.9, 3.5};
	float lens_y[5] = {4, 3.5, 1.75, 2.6, 3.4};
	float lens_mass[5] = {1, 2, 3, 4, 5};
	/* end test data */
	
	cell ** cells = malloc(sizeof(cell) * max_index);
	curr = (cell *)malloc(sizeof(cell));
	curr->index = 0;
	curr->top_left[0] = 0;
	curr->top_left[1] = pixel_y;
	curr->bottom_right[0] = pixel_x;
	curr->bottom_right[1] = 0;
	
	root = curr;
	cells[0] = root;
	
	for (i = 0; i < nobjects; ++i) {
		dx = lens_x[i];
		dy = lens_y[i];
		lmass = lens_mass[i];
		printf("Current lens: %f %f\n", dx, dy);
		
		// if lens is contained within the dimensions of the current cell
		while (dx > curr->top_left[0] && dx < curr->bottom_right[0] &&
			   dy > curr->bottom_right[1] && dy < curr->top_left[1]) {
			
			printf("Cell is within range\n");
			
			// If there are no lens and subcells in the cell, set curr to be the lens, otherwise move to the child cell. 
			if (!curr->clens.pos_x && !curr->subcells[0]) {
				printf("No lens in current cell\n");
				curr->clens.pos_x = dx;
				curr->clens.pos_y = dy;
				curr->clens.mass = lmass;
				
				// TO DO: set other cell parameters too
				break;
			
			} else {
				// More than one lens being placed in a cell
				// If there are subcells, find which of the subcells cells the lens belongs in, otherwise divide cell.
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
					
					do {
						lens1_index = 0;
						lens2_index = 0;
						for (j = 0; j < ncells; ++j) {
							temp = (cell *)malloc(sizeof(cell));
							
							// calculate the index of the array in which the lens will be stored
							temp->index = ncells * curr->index + (j+1);
							curr->subcells[j] = temp->index;
							
							// work out the dimensions of the subcell
							temp->top_left[0] = (curr->bottom_right[0] - curr->top_left[0])/2 * (j % 2) + curr->top_left[0];
							temp->top_left[1] = j > 1 ? (curr->top_left[1] - curr->bottom_right[1])/2 + curr->bottom_right[1] : curr->top_left[1];
							temp->bottom_right[0] = (j%2) == 0 ? (curr->bottom_right[0] - curr->top_left[0])/2 + curr->top_left[0] : curr->bottom_right[0];
							temp->bottom_right[1] = j < 2 ? (curr->top_left[1] - curr->bottom_right[1])/2 + curr->bottom_right[1] : curr->bottom_right[1];
							
							cells[temp->index] = temp;
						}
						
						// put the lens currently in the cell to be divided and the new lens in the subcells. 
						for (j = 0; j < ncells; ++j) {
							temp = cells[curr->subcells[j]];
							
							if (dx > temp->top_left[0] && dx < temp->bottom_right[0] &&
								dy > temp->bottom_right[1] && dy < temp->top_left[1]) {
								temp->clens.pos_x = dx;
								temp->clens.pos_y = dy;
								temp->clens.mass = lmass;
								lens1_index = temp->index;
							}
							
							if (curr->clens.pos_x > temp->top_left[0] && curr->clens.pos_x < temp->bottom_right[0] &&
								curr->clens.pos_y > temp->bottom_right[1] && curr->clens.pos_y < temp->top_left[1]) {
								temp->clens = curr->clens;
								curr->clens.pos_x = 0;
								curr->clens.pos_y = 0;
								curr->clens.mass = 0;
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
	
	printf("\nIndex of cell and its dimensions. i.e. top left corner and bottom right corner x, y co-ordinates\n");
	for (i = 0; i < max_index; ++i) {
		if(cells[i]) {
		printf("%d %f %f %f %f\n", i, cells[i]->top_left[0], cells[i]->top_left[1], cells[i]->bottom_right[0], cells[i]->bottom_right[1]);
		printf("lens: %f %f %f\n", cells[i]->clens.pos_x, cells[i]->clens.pos_y, cells[i]->clens.mass);
		}
	}
	
	free(curr);
	// free cell array
	return 0;
}

void calculate_r(cell ** cells) {
	printf("calculate the centre of mass and add up masses");
}
