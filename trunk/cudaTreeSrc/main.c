#include "global.h"
#include "constants.h"
#include "util.h"
#include "tree_struct.h"

#define pixel_x (512)
#define pixel_y (512)

//float4 *lenses;
float *lens_x;
float *lens_y;
float *lens_mass;
size_t nobjects;
cell *root; 	/* global pointer to root cell so we can get to it from anywhere */
const float delta = 0.5;		/*variable to determine accuracy*/

int main(int argc, const char *argv[]){
  int i;
  root = NULL;
  //cell *temp;
  lens *temp_lens;
  // Load relevant settings and data
  if (argc < 2) error("Requires argument with lens positions and optional mass");
  setup_constants();
  read_lenses(argv[1]);
  make_root(&root);
  printf("Values from root node: width: %f. height: %d, index: %d, top-left: %f %f, bot_right: %f %f\n", root->width, root->height, root->index, root->top_left[0], root->top_left[1], root->bot_right[0], root->bot_right[1]);
  for(i=0; i < nobjects; i++){
    //temp = (cell *)salloc(sizeof(cell));                /*allocating memory to temp cell*/
    temp_lens = (lens *)salloc(sizeof(lens));                /*allocating memory to temp lens*/
    temp_lens->x = lens_x[i];                          /*assigning x-coord to cell_lens x-coord*/
    temp_lens->y = lens_y[i];                          /*assigning y-coord to cell_lens y-coord*/
    temp_lens->m = lens_mass[i];                          /*assigning mass*/
    //printf("before we call the make_tree method\n");
    make_tree(&root, temp_lens);
    free(temp_lens);
  }
  printf("The total center of mass for all the lenses is: %.4f\n", root->total_mass);
  printf("The total value for the numerator for x is: %.4f\n", root->center_mass_x);
  printf("The center of mass in x direction is: %.4f\n", root->center_mass_x/root->total_mass);
  printf("The total value for the numerator for y is: %.4f\n", root->center_mass_y);
  printf("The center of mass in y direction is: %.4f\n", root->center_mass_y/root->total_mass);
  printf("The x and y value for the top left corner is: x:%.4f, y: %.4f\n", root->top_left[0], root->top_left[1]);
  printf("The x and y value for the bottom right corner is: x:%.4f, y: %.4f\n", root->bot_right[0], root->bot_right[1]);

  printf("printing values from tree: \n");
  print_tree(&root);
  
  free_tree(&root);
  return 0;
}
