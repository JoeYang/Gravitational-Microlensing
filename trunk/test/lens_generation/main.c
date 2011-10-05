#include <stdio.h>
#include "lens_gen.c"

#define pixel	20
#define rpp		20

int main(){
	int size = pixel*pixel*rpp;
	float *lens_gen_x = (float*)malloc(size*sizeof(float));
	float *lens_gen_y = (float*)malloc(size*sizeof(float));
	
	lens_generate(lens_gen_x, lens_gen_y, rpp, pixel);	

	int it;
	for(it = 0; it< size; it++){
		printf("%f %f\n", lens_gen_x[it], lens_gen_y[it]);
	}
	
	return 0;
}