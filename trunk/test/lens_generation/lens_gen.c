#include "constants.h"
#include <math.h>

void lens_generate(float *new_lens_x, float *new_lens_y, int rpp, int pixel){
	
	setup_constants();
	float increment_x = (image_scale_x * 2) / pixel;
  	float increment_y = (image_scale_y * 2) / pixel;
	
	int size = pixel*pixel*rpp;
	int i = 0;
	int it;	
	srand(time(NULL));
	
	for(it = 0; it < rpp; ++it){
		float x, y;
		for(y = -image_scale_y; y < image_scale_y; y += increment_y) {
      		for(x = -image_scale_x; x < image_scale_x; x += increment_x) {
		        // Noise is uniformly distributed -- i.e. it's not really noise
		        float noise_x = increment_x*rand() / RAND_MAX;
		        float noise_y = increment_x*rand() / RAND_MAX;
		        new_lens_x[i] = x + noise_x;
		        new_lens_y[i] = y + noise_y;
  			   	i++;

            }
    	} 
	}
}
