#include "tree_code.h"

int main(void) {
	cell * curr, * root, * temp;
	int i, j, index;
	float dx, dy;
	int lmass;
	int max_index = 20; // TO DO: need to work out max number of indexes that can potentially be used
	cell * cells[max_index];
	int ncells = 4; // number of cells to split a cell into when two lenses are in the same cell
	int nobjects = 2;
	
	int pixel_x = 5, pixel_y = 5; // dimensions of the lensing plane
	float lens_x[2] = {2, 3.5};
	float lens_y[2] = {4, 3.5};
	float lens_mass[2] = {1, 2};
	
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
			
			// If there are no lens and children cells in the cell, set curr to be the lens, otherwise move to the child cell. 
			if (!curr->clens.pos_x && !curr->children[0]) {
				printf("No lens in current cell\n");
				curr->clens.pos_x = dx;
				curr->clens.pos_y = dy;
				curr->clens.mass = lmass;
				
				// TO DO: set other cell parameters too
				break;
			
			} else {
				// More than one lens being placed in a cell
				// If there are children, find which of the children cells the lens belongs in, otherwise divide cell.
				if (curr->children[0]) {
					printf("Found children\n");
					
					for (j = 0; j < ncells; ++j) {
						temp = cells[curr->children[j]];
						if (dx > temp->top_left[0] && dx < temp->bottom_right[0] &&
							dy > temp->bottom_right[1] && dy < temp->top_left[1]) {
							
							curr = temp;
							break;
						}
					}
				} else {
					printf("No children. Need to create new cells\n");
					
					for (j = 0; j < ncells; ++j) {
						temp = (cell *)malloc(sizeof(cell));
						
						// calculate the index of the array in which the lens will be stored
						temp->index = 4 * curr->index + (j+1);
						curr->children[j] = temp->index;
						
						// work out the dimensions of the child cell
						temp->top_left[0] = (curr->bottom_right[0] - curr->top_left[0])/2 * (j % 2) + curr->top_left[0];
						temp->top_left[1] = j > 1 ? (curr->top_left[1] - curr->bottom_right[1])/2 + curr->bottom_right[1] : curr->top_left[1];
						temp->bottom_right[0] = (j%2) == 0 ? (curr->bottom_right[0] - curr->top_left[0])/2 + curr->top_left[0] : curr->bottom_right[0];
						temp->bottom_right[1] = j < 2 ? (curr->top_left[1] - curr->bottom_right[1])/2 + curr->bottom_right[1] : curr->bottom_right[1];
						
						cells[temp->index] = temp;
					}
					
					// put the lens in the divided cell and the new lens in the child cells. 
					for (j = 0; j < ncells; ++j) {
						temp = cells[curr->children[j]];
						
						if (dx > temp->top_left[0] && dx < temp->bottom_right[0] &&
							dy > curr->bottom_right[1] && dy < temp->top_left[1]) {
							
							temp->clens.pos_x = dx;
							temp->clens.pos_y = dy;
							temp->clens.mass = lmass;
						}
						
						if (curr->top_left[0] > temp->top_left[0] && curr->bottom_right[0] < temp->bottom_right[0] &&
							curr->bottom_right[1] > temp->bottom_right[1] && curr->top_left[1] < temp->top_left[1]) {
							
							temp->clens = curr->clens;
							curr->clens.pos_x = 0;
							curr->clens.pos_y = 0;
							curr->clens.mass = 0;
						}
						
						// TO DO: Need to consider if the two lenses are still in the same cell after the cells are divided. Recursion is looking more and more like a better option. :(
					}
					
					break;
					
				}
			}
		}
		
		// when lens comes here it either lies out of the grid or it was successfully assigned
		// TO DO: Deal with lenses that do not fit in the grid
		
		curr = root;
	}
	
	printf("\nIndex of cell and its dimensions. i.e. top left corner and bottom left corner x, y co-ordinates\n");
	for (i = 0; i < 5; ++i) {
		printf("%d %f %f %f %f\n", i, cells[i]->top_left[0], cells[i]->top_left[1], cells[i]->bottom_right[0], cells[i]->bottom_right[1]);
	}
	
	free(curr);
	// free all cells
	return 0;
}
