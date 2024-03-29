#ifndef TREE_STRUCT_H
#define TREE_STRUCT_H

#define DIM 2
#define DESC (1 << DIM) /*calculate the number of decendents a cell has*/
/*basic struct to hold data of lens and pointer
 * that is used by both cells and */
typedef struct lens{
  float x;		/*x position of the lens*/
  float y;		/*y position of the lens*/
  float m;		/*mass of the lens*/
  //struct lens *next;	/*pointer to the next node (may not need this)*/
}lens;

/* struct of a cell in the tree that has
 * pointers to other cells or lens in the tree*/
typedef struct cell{
  int index;
  float width;
  float length;
  int height;
  float total_mass;	/*total mass, but can be used as center of mass denomerator: sum(lens_mass[i])*/
  float center_mass_x;	/*center of mass numerator of x-pos: sum(lens_mass*lens_x[i])*/
  float center_mass_y;	/*center of mass numerator of 'y'-pos: sum(mass*x-pos)*/
  float quad_pole;
  float multipole;
  float top_left[2];
  float bot_right[2];
  lens *lenses[DESC]; 	/* an array of the lenses of the cell*/
  struct cell *desc[DESC]; /* an array of descendents of the cell*/
}cell;

void make_root(cell **root);
void make_tree(cell ** tree, lens *lens);
void print_tree(cell ** tree);
//void get_lens_count(cell ** tree, float ray_x, float ray_y);
//int get_lens_count(cell ** tree, float delta, float ray_x, float ray_y);
void get_lens_count(cell ** tree, float delta, float ray_x, float ray_y, int * lens_count);
void free_tree(cell ** tree);

#endif
