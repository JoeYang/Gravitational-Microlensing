#include "global.h"

#ifndef CONSTANTS_H_
#define CONSTANTS_H_

extern float kappa_star, kappa_c, gamma_, kappa, source_scale;
extern float image_scale_fudge, image_scale_x, image_scale_y;
extern float lens_scale_fudge_1, lens_scale_fudge_2, lens_rad_x, lens_rad_y, lens_rad;
extern float lens_count;

void setup_constants(void);

#endif
