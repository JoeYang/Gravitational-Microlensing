#ifndef TREECODE_H_
#define TREECODE_H_

#include <stdio.h>
#include <stdlib.h>

/* data structure to store a lens' position and mass */
typedef struct lens {
	float pos_x;
	float pos_y;
	float mass;
} lens;

/* data structure to store information about a cell */
typedef struct cell {
	int index; // position where cell is stored in array, 'cells'
	float top_left[2], bottom_right[2]; // x, y co-ordinates that define the dimensions of the cell
	float centre_of_mass;
	float mass; // total mass of all lenses in the cell
	float quadrupole_moment;
	float higher_multipole_moment; 
	int children[4]; // indexes of the subcells
	lens clens; // lens inside the cell. If the lens has subcells this doesn't exist.
} cell;

#endif