#include <stdio.h>
#include "lens_gen.c"

#define size	20

int main(){
	
	float *lens_gen_x = (float*)malloc(size*sizeof(float));
	float *lens_gen_y = (float*)malloc(size*sizeof(float));
	
	lens_generate(lens_gen_x, lens_gen_y, size);	

	int it;
	for(it = 0; it< size; it++){
		printf("%f %f\n", lens_gen_x[it], lens_gen_y[it]);
	}
	
	return 0;
}