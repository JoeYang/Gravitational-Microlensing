#include "constants.h"
#include <math.h>

void lens_generate(float *new_lens_x, float *new_lens_y, int nobjects){
	
	setup_constants();
		
	int it;	
	srand(time(NULL));
	
	for(it = 0; it < nobjects; ++it){
		float x, y;
        x = lens_rad_x*rand() / RAND_MAX;
        y = lens_rad_y*rand() / RAND_MAX;
        new_lens_x[it] = x;
        new_lens_y[it] = y; 
	}
}
