#ifndef TREECODE_H_
#define TREECODE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* data structure to store a lens' position and mass */
typedef struct lens lens;
typedef struct cell cell;

struct lens {
	float pos_x;
	float pos_y;
	float mass;
};

/* data structure to store information about a cell */
struct cell {
	int index; // position where cell is stored in array, 'cells'
	float top_left[2], bottom_right[2]; // x, y co-ordinates that define the dimensions of the cell
	float cm_x; // centre of mass, x
	float cm_y; // centre of mass, y
	float mass; // total mass the cell's lenses
	cell *subcells[4]; // indexes of the subcells
	lens clens; // lens inside the cell. If the lens has subcells this doesn't exist.
	/* float quadrupole_moment; */
	/* float higher_multipole_moment; */
};

#endif