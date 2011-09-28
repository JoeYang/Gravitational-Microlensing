#ifndef TREECODE_H_
#define TREECODE_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
	float cm_x; // centre of mass, x
	float cm_y; // centre of mass, y
	float mass; // total mass the cell's lenses
	struct cell *subcells[4]; // indexes of the subcells
	struct lens clens; // lens inside the cell. If the lens has subcells this doesn't exist.
	/* float quadrupole_moment; */
	/* float higher_multipole_moment; */
} cell;

//void setup_root(float, float, float, float, int);
void build_tree(int, cell *);
void calculate_cm(cell *, int);
//void *salloc(size_t);
void print_tree(cell *, int);
void free_tree(cell *, int);
void remove_empty_cells(cell *, int);
void get_included_bodies(cell *, float, float, float);

#endif
