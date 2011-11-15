#ifndef TREECODE_KERNEL_H
#define TREECODE_KERNEL_H

/* structure to hold the variables so we can pass them to the GPU*/
typedef struct vars {
  unsigned int rpp;
  float kappa_c, gamma_, source_scale;
  float image_scale_x, image_scale_y;
  float increment_x, increment_y;
} vars;

/*void cudaHandleError( cudaError_t err, const char *file, int line){
  if( err != cudaSuccess ){
    printf("%s in %s at line %d\n", cudaGetErrorString(err), file, line);
    exit( EXIT_FAILURE );
  }
}*/

#define HANDLE_ERROR( err ) (cudaHandleError(err, __FILE__, __LINE__))

//void init_var(vars *var);
//int highest(unsigned int *results, unsigned int size);
__global__ void glensing(const float4 * lenses, const size_t nobjects, unsigned int* results, const vars* v);

#endif
