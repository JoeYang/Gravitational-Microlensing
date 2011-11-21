#include "tree_struct.h"
#include "util.h"
#include "constants.h"
#include "global.h"

#define QUAD_IDX (4)

/* method to create root of tree,
 * easier to create root separately
 * so we dont have to worry about null variables
 */
void make_root(cell **root){
  /* adding the values of the top left and bottom right
   * of the node plane to the head of the tree
   */
  int i;
  (*root) = (cell *)salloc(sizeof(cell));
  //(*root)->width = 1.0;
  (*root)->length = 1.0;
  (*root)->index = 0;
  (*root)->height = 0;
  (*root)->top_left[0] = -1.0*(float)(lens_rad);
  (*root)->top_left[1] = (float)(lens_rad);
  (*root)->bot_right[0] = (float)(lens_rad);
  (*root)->bot_right[1] = -1.0*(float)(lens_rad);
  (*root)->width = (*root)->bot_right[0] - (*root)->top_left[0];
  for(i=0; i < QUAD_IDX; i++){
    (*root)->lenses[i] = NULL;
    (*root)->desc[i] = NULL;
  }
}

/* method for getting the number of lenses to include
*  if the ratio is greater than one we resolve it to one
*
*/
void get_lens_count(cell ** tree, float delta, float ray_x, float ray_y, int * lens_count){
  //static int val = 0;                 /* static variable to contain the number of lenses to include */
  int i, j, has_cells=0, has_lenses=0;

  /* finding the values of the current cell to decide if it is to be included in the calculation*/
  float cell_x = (float)(*tree)->center_mass_x/(float)(*tree)->total_mass;
  float cell_y = (float)(*tree)->center_mass_y/(float)(*tree)->total_mass;
  float dist = sqrt(pow(ray_x - cell_x,2) + pow(ray_y - cell_y, 2));
  float ratio = ((*tree)->width)/dist;

  // checking if the current cell has any cells or lenses
  for(j=0; j<QUAD_IDX; j++){
    if((*tree)->lenses[j]) has_lenses++;
    if((*tree)->desc[j]) has_cells++;
  }

  if(delta == 0.0 && has_lenses>0){
    (*lens_count) += has_lenses;
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->desc[i]) get_lens_count(&((*tree)->desc[i]), delta, ray_x, ray_y, lens_count);
    }
  }
  else if(has_lenses==1 && has_cells==0){
    (*lens_count)++;
    return;
  }
  else if(has_lenses>0 && has_cells==0 && ratio>delta){
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->lenses[i]){
        (*lens_count)++;
      }
    }
    return;
  }
  else if(has_lenses>0 && has_cells>0 && ratio>delta){
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->lenses[i]){
        (*lens_count)++;
      }
      if((*tree)->desc[i]) get_lens_count(&((*tree)->desc[i]), delta, ray_x, ray_y, lens_count);
    }
    return;
  }
  else if(ratio<delta){
    (*lens_count)++;
    return;
  }
  else{
    int i;
    for(i=0; i<QUAD_IDX; i++){
      if((*tree)->desc[i]) get_lens_count(&((*tree)->desc[i]), delta, ray_x, ray_y, lens_count);
    }
  }
  //return (*lens_count);
}

/* 
 * Method to return which quadrant the lens is in:
 * 0 - for top-right
 * 1 - for top-left
 * 2 - for bottom-left
 * 3 - for bottom-right
 */
int get_lens_pos(cell *tree, lens *lens){
  //this checks to see if the lens is in the top right quadrant
  if(lens->x >= (tree->top_left[0] + tree->bot_right[0])/2.0 && lens->x <= tree->bot_right[0]
          && lens->y >= (tree->top_left[1] + tree->bot_right[1])/2.0  && lens->y <= tree->top_left[1]){
    return 0;
  }
  //this check is to see if the lens is in the top left quadrant
  else if(lens->x >= tree->top_left[0] && lens->x <= (tree->top_left[0] + tree->bot_right[0])/2.0
     && lens->y >= (tree->top_left[1] + tree->bot_right[1])/2.0 && lens->y <= tree->top_left[1]){
    return 1;
  }
  //check if the lens is in the bottom left quadrant
  else if(lens->x >= tree->top_left[0] && lens->x <= (tree->top_left[0] + tree->bot_right[0])/2.0
           && lens->y >= tree->bot_right[1] && lens->y <= (tree->top_left[1] + tree->bot_right[1])/2.0){
    return 2;
  }
  //check if the lens is in the bottom right quadrant (maybe change this later to an else)
  else if(lens->x >= (tree->top_left[0] + tree->bot_right[0])/2.0 && lens->x <= tree->bot_right[0]
          && lens->y >= tree->bot_right[1] && lens->y <= (tree->top_left[1] + tree->bot_right[1])/2.0){
    return 3;
  }
  else{
   // printf("The x and y values from the lens: x: %.8f, y: %.8f\n", lens->x, lens->y);
   // printf("The x and y values from the top left of tree: x: %.8f, y: %.8f\n", tree->top_left[0], tree->top_left[1]);
   // printf("The x and y values from the bottom right of tree: x: %.8f, y: %.8f\n", tree->bot_right[0], tree->bot_right[1]);
    error("could not find the lens position\n");
  }

}

/*
 * Method to create a cell when needed by make_tree
 **/
cell * create_cell(cell ** tree, lens *nlens, lens *olens, int pos){
  int i, count = 0, recip;          /* The reciprocal of the height of the cell */
  cell *c = (cell *)salloc(sizeof(cell));
  c->total_mass = olens->m + nlens->m;      /* total mass of the cell*/
  c->center_mass_x = (olens->m * olens->x) + (nlens->m * nlens->x); /* numerator for calculation of center of mass in x */
  c->center_mass_y = (olens->m * olens->y) + (nlens->m * nlens->y); /* numerator for calculation of center of mass in y */
  //c->width = (*tree)->width/2.0;                                    /* width of cell in the graph*/
  c->length = (*tree)->length/2.0;                                    /* width of cell in the graph*/
  c->index = QUAD_IDX*((*tree)->index) + pos + 1;               /* way to calculate index of node in tree*/
  recip = (int)(1.0/(float)(c->length));
  while(recip > 1){
    recip /= 2;
    count++;
  }
  c->height = count;
  for(i=0; i < QUAD_IDX; i++){          /* setting all lenses  and descendents of the current node to NULL*/
    c->lenses[i] = NULL;
    c->desc[i] = NULL;
  }
  /*
  * the quadrants are layed out as such:
  *  _________________________
  * |            |            |
  * | Quadrant 2 | Quadrant 1 |
  * |____________|____________|
  * |            |            |
  * | Quadrant 3 | Quadrant 4 |
  * |____________|____________|
  */
  /* finding the correct position for the cell*/
  if(pos == 0){               /* first quadrant */
    c->top_left[0] = ((*tree)->top_left[0] + (*tree)->bot_right[0])/2.0;
    c->top_left[1] = (*tree)->top_left[1];
    c->bot_right[0] = (*tree)->bot_right[0];
    c->bot_right[1] = ((*tree)->bot_right[1] + (*tree)->top_left[1])/2.0;
  }
  else if(pos == 1){          /* second quadrant */
    c->top_left[0] = (*tree)->top_left[0];
    c->top_left[1] = (*tree)->top_left[1];
    c->bot_right[0] = ((*tree)->bot_right[0] + (*tree)->top_left[0])/2.0;
    c->bot_right[1] = ((*tree)->bot_right[1] + (*tree)->top_left[1])/2.0;
  }
  else if(pos == 2){          /* third quadrant */
    c->top_left[0] = (*tree)->top_left[0];
    c->top_left[1] = ((*tree)->top_left[1] + (*tree)->bot_right[1])/2.0;
    c->bot_right[0] = ((*tree)->bot_right[0] + (*tree)->top_left[0])/2.0;
    c->bot_right[1] = (*tree)->bot_right[1];
  }
  else if(pos == 3){          /* fourth quadrant */
    c->top_left[0] = ((*tree)->top_left[0] + (*tree)->bot_right[0])/2.0;
    c->top_left[1] = ((*tree)->top_left[1] + (*tree)->bot_right[1])/2.0;
    c->bot_right[0] = (*tree)->bot_right[0];
    c->bot_right[1] = (*tree)->bot_right[1];
  }
  c->width = c->bot_right[0] - c->top_left[0];

  return c;
}

/*
 * Method to create a new cell from the two lenses that
 * we are trying to put into the same cell
 */
void add_lens(cell ** tree, lens **olens, lens **nlens, int pos){
  int n_lens_pos = get_lens_pos(*tree, *nlens);                   /* The quadrant of the lens we want to insert*/
  int o_lens_pos = get_lens_pos(*tree, *olens); /* The quadrant of the lens that already exists*/

  if(!(n_lens_pos == o_lens_pos)){
    (*tree)->lenses[o_lens_pos] = *olens;
    (*tree)->lenses[n_lens_pos] = *nlens;
  }
  else{
    /* if we get here the positions are the same for the lenses, but maybe different from the previous position obtained*/
    int new_pos = n_lens_pos;
    cell *tmp = create_cell(tree, *olens, *nlens, new_pos);
    add_lens(&tmp, olens, nlens, new_pos);
    (*tree)->desc[new_pos] = tmp;
    (*tree)->lenses[new_pos] = NULL;
  }
}

/*method to create Barnes-Hut tree*/
void make_tree(cell ** tree, lens *lens){
  int i, lens_pos;
  struct lens *temp = (struct lens *)salloc(sizeof(struct lens));                /*allocating memory to temp lens*/
  /* adding the properties of *lens to *temp */
  temp->x = lens->x;
  temp->y = lens->y;
  temp->m = lens->m;
  lens_pos = get_lens_pos(*tree, lens);
  (*tree)->total_mass += lens->m;
  (*tree)->center_mass_x += lens->m*lens->x;
  (*tree)->center_mass_y += lens->m*lens->y;

	if(!((*tree)->desc[0]) && !((*tree)->desc[1]) && !((*tree)->desc[2]) && !((*tree)->desc[3])){
    
    if(!((*tree)->lenses[lens_pos])){
      (*tree)->lenses[lens_pos] = temp;
    }
    else if(((*tree)->lenses[lens_pos]->x == temp->x) && ((*tree)->lenses[lens_pos]->y == temp->y)){
      (*tree)->lenses[lens_pos]->m += temp->m;
    }
    else{
      /* we have a lens here and so we need to create a new descendent node*/
      add_lens(tree, &((*tree)->lenses[lens_pos]), &temp, lens_pos);
      (*tree)->lenses[lens_pos] = NULL;
    }
      
  }
  else{     /* here if there are already cells */

    /* if there is a cell in this position*/
    if(!((*tree)->desc[lens_pos]) && !((*tree)->lenses[lens_pos])){
      (*tree)->lenses[lens_pos] = temp;
    }
    else if((*tree)->desc[lens_pos] && !((*tree)->lenses[lens_pos])){
      free(temp);
      make_tree(&((*tree)->desc[lens_pos]), lens);
    }
    else if(!((*tree)->desc[lens_pos]) && (*tree)->lenses[lens_pos]){
      if(((*tree)->lenses[lens_pos]->x == temp->x) && ((*tree)->lenses[lens_pos]->y == temp->y)){
        (*tree)->lenses[lens_pos]->m += temp->m;
        free(temp);
      }
      else{
        add_lens(tree,&((*tree)->lenses[lens_pos]), &temp, lens_pos);
        (*tree)->lenses[lens_pos] = NULL;
      }
    }
    else{
      error("Error occurred while constructing tree");
    }
  }
}

/* method to print tree*/
void print_tree(cell ** tree){
  int i;
  printf("\n");
  printf("Values of Cell:\n");
  printf("Index: %d\n", (*tree)->index);
  printf("Height: %d\n", (*tree)->height);
  printf("Width: %.8f\n", (*tree)->width);
  printf("center of mass x: %.8f\n", ((*tree)->center_mass_x)/((*tree)->total_mass));
  printf("center of mass y: %.8f\n", ((*tree)->center_mass_y)/((*tree)->total_mass));
  printf("values for any lenses the cell has: \n");
  for(i=0; i<QUAD_IDX; i++){
    if((*tree)->lenses[i]){
      printf("X, Y and mass for lens in position: %d\n", i);
      printf("lens x-coord: %.8f\n", (*tree)->lenses[i]->x);
      printf("lens y-coord: %.8f\n", (*tree)->lenses[i]->y);
    }
  }
  printf("printing the values : \n");
  for(i=0; i<QUAD_IDX; i++){
    if((*tree)->desc[i]) print_tree(&((*tree)->desc[i]));
  }
}

/* Method to free memory used by tree*/
void free_tree(cell ** tree){

  if(!((*tree)->desc[0]) && !((*tree)->desc[1]) && !((*tree)->desc[2]) && !((*tree)->desc[3])){
    int i;
    for(i=0; i< QUAD_IDX; i++){
      if((*tree)->lenses[i]) free((*tree)->lenses[i]);
    }
    free(*tree);
  }
  else{
    int i;
    for(i=0; i< QUAD_IDX; i++){
      if((*tree)->lenses[i]) free((*tree)->lenses[i]);
      if((*tree)->desc[i]) free_tree(&((*tree)->desc[i]));
    }
    free(*tree);
  }
}
