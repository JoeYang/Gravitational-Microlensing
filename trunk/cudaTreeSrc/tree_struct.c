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
  (*root)->width = 1.0;
  (*root)->index = 0;
  (*root)->height = 0;
  (*root)->top_left[0] = -1.0*(float)(lens_rad);
  (*root)->top_left[1] = (float)(lens_rad);
  (*root)->bot_right[0] = (float)(lens_rad);
  (*root)->bot_right[1] = -1.0*(float)(lens_rad);
  for(i=0; i < QUAD_IDX; i++){
    (*root)->lenses[i] = NULL;
    (*root)->desc[i] = NULL;
  }
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
  else if(lens->x >= tree->top_left[0] && lens->x < (tree->top_left[0] + tree->bot_right[0])/2.0
     && lens->y > (tree->top_left[1] + tree->bot_right[1])/2.0 && lens->y <= tree->top_left[1]){
    return 1;
  }
  //check if the lens is in the bottom left quadrant
  else if(lens->x >= tree->top_left[0] && lens->x < (tree->top_left[0] + tree->bot_right[0])/2.0
           && lens->y >= tree->bot_right[1] && lens->y < (tree->top_left[1] + tree->bot_right[1])/2.0){
    return 2;
  }
  //check if the lens is in the bottom right quadrant (maybe change this later to an else)
  else if(lens->x >= (tree->top_left[0] + tree->bot_right[0])/2.0 && lens->x <= tree->bot_right[0]
          && lens->y >= tree->bot_right[1] && lens->y <= (tree->top_left[1] + tree->bot_right[1])/2.0){
    return 3;
  }
  else{
 //   printf("The x and y values from the lens: x: %.8f, y: %.8f\n", lens->x, lens->y);
 //   printf("The x and y values from the top left of tree: x: %.8f, y: %.8f\n", tree->top_left[0], tree->top_left[1]);
 //   printf("The x and y values from the bottom right of tree: x: %.8f, y: %.8f\n", tree->bot_right[0], tree->bot_right[1]);
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
  c->width = (*tree)->width/2.0;                                    /* width of cell in the graph*/
  c->index = QUAD_IDX*((*tree)->index) + pos + 1;               /* way to calculate index of node in tree*/
  recip = (int)(1.0/(float)(c->width));
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
  //printf("==================================================================================================\n");
  int i, lens_pos;
  struct lens *temp = (struct lens *)salloc(sizeof(struct lens));                /*allocating memory to temp lens*/
  temp->x = lens->x;
  temp->y = lens->y;
  temp->m = lens->m;
  lens_pos = get_lens_pos(*tree, lens);
  //printf("the position of the current lens: %d\n", lens_pos);
  //printf("The values of the current lens: x: %f, y: %f\n", lens->x, lens->y);
  (*tree)->total_mass += lens->m;
  (*tree)->center_mass_x += lens->m*lens->x;
  (*tree)->center_mass_y += lens->m*lens->y;
  //printf("value of the temp node: x: %f, y: %f, m: %f\n", temp->x, temp->y, temp->m);
  //printf("The positions of the tree node: top_left x: %f, top_left y: %f\n", (*tree)->top_left[0], (*tree)->top_left[1]);
  //printf("The positions of the tree node: bot_right x: %f, bot_right y: %f\n", (*tree)->bot_right[0], (*tree)->bot_right[1]);

	if(!((*tree)->desc[0]) && !((*tree)->desc[1]) && !((*tree)->desc[2]) && !((*tree)->desc[3])){
    //printf("current tree node has no descendents\n");
    
    if(!((*tree)->lenses[lens_pos])){
      (*tree)->lenses[lens_pos] = temp;
    }
    else{
      /* we have a lens here and so we need to create a new descendent node*/
      //printf("The values of the lens in the tree: x: %f, y: %f, m: %f\n", (*tree)->lenses[lens_pos]->x, (*tree)->lenses[lens_pos]->y, (*tree)->lenses[lens_pos]->m);
      add_lens(tree, &((*tree)->lenses[lens_pos]), &temp, lens_pos);
      (*tree)->lenses[lens_pos] = NULL;
    }
      
  }
  else{     /* here if there are already cells */
    //printf("There are already cells in the tree\n");

    /* if there is a cell in this position*/
    if(!((*tree)->desc[lens_pos]) && !((*tree)->lenses[lens_pos])){
      (*tree)->lenses[lens_pos] = temp;
    }
    else if((*tree)->desc[lens_pos] && !((*tree)->lenses[lens_pos])){
      free(temp);
      //printf("The top_left of the descendent: %d x: %.8f, y: %.8f\n", lens_pos, (*tree)->desc[lens_pos]->top_left[0], (*tree)->desc[lens_pos]->top_left[1]);
     // printf("The bottom_right of the descendent: %d x: %.8f, y: %.8f\n", lens_pos, (*tree)->desc[lens_pos]->bot_right[0], (*tree)->desc[lens_pos]->bot_right[1]);
      make_tree(&((*tree)->desc[lens_pos]), lens);
    }
    else if(!((*tree)->desc[lens_pos]) && (*tree)->lenses[lens_pos]){
      //printf("The values of the lens in the tree: x: %f, y: %f, m: %f\n", (*tree)->lenses[lens_pos]->x, (*tree)->lenses[lens_pos]->y, (*tree)->lenses[lens_pos]->m);
      add_lens(tree,&((*tree)->lenses[lens_pos]), &temp, lens_pos);
      (*tree)->lenses[lens_pos] = NULL;
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
  printf("values for any lenses the cell has: \n");
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
