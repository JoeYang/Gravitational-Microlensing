 #include "constants.h"

float kappa_star, kappa_c, gamma_, kappa, source_scale;
float image_scale_fudge, image_scale_x, image_scale_y;
float lens_scale_fudge_1, lens_scale_fudge_2, lens_rad_x, lens_rad_y, lens_rad;
float lens_count;

vars* setup_constants() {
	/* Two parameters for the lensing system: convergence (kappa) and shear (gamma).
	   Convergence can be broken into a smooth component (kappa_c) and a clumpy stellar
	   component (kappa_star), such that kappa = kappa_star + kappa_c */
	kappa_star = 0.394;           // convergence in stars (sometimes written as sigma instead of kappa)
	kappa_c = 0.;                 // convergence in smooth matter
	gamma_ = 0.395;                // shear
	kappa = kappa_star + kappa_c; // total convergence
	/* Source plane (where the light rays are collected) is a square with side length source_scale,
	   measured in Einstein Radii. */
	source_scale = 24.;       // size of the source plane in Einstein Radii
	/* Size of the image plane (from which the light rays are shot) depends on the desired size of the 
	   source plane and the lensing parameters. Note that the image plane is rectangular -- 2d ray 
	   positions should be randomly generated in the ranges [-image_scale_x, image_scale_x] and
	   [-image_scale_y, image_scale_y].
	   Fudge factor values were determined empirically -- leaving them set to the following values 
	   will allow direct comparison of output with single core tree code */
	image_scale_fudge = 0.1;  // image plane scale fudge factor
	image_scale_x = (0.5 + 2.0*image_scale_fudge) * source_scale / (1.0 - kappa - gamma_);
	image_scale_y = (0.5 + 2.0*image_scale_fudge) * source_scale / (1.0 - kappa + gamma_);
	/* Radius of the lens plane (where the microlenses are distributed) also depends on the desired size
	   of the source plane and the lensing parameters. We use a circular lens plane of radius lens_rad.
	   Some notes apply about the fudge factors. */
	lens_scale_fudge_1 = 5.0; // lens plane scale fudge factor 1
	lens_scale_fudge_2 = 2.0; // lens plane scale fudge factor 2
	lens_rad_x = (0.5*source_scale + lens_scale_fudge_1) / (1.0 - kappa + gamma_);
	lens_rad_y = (0.5*source_scale + lens_scale_fudge_1) / (1.0 - kappa - gamma_);
	lens_rad = sqrt(lens_rad_x*lens_rad_x + lens_rad_y*lens_rad_y) + lens_scale_fudge_2;
	/* lens_count is the required number of microlenses. These should be uniformly randomly distributed 
	   within a circle of radius lens_rad. */  
	lens_count = (size_t)(kappa_star * lens_rad*lens_rad / 1 + 0.5);
	// TODO: Eventually this should be the avg. lens mass
	// 1 is the lens mass
 	vars *var = (vars*)malloc(sizeof(vars));
 	init_var(var);
	return var;
}

void init_var(vars *var){
	var->kappa_star = kappa_star;
	var->kappa_c = kappa_c;
	var->gamma_ = gamma_;
	var->kappa = kappa; 
	var->source_scale = source_scale;
	var->image_scale_fudge = image_scale_fudge;
	var->image_scale_x = image_scale_x; 
	var->image_scale_y = image_scale_y;
	var->lens_scale_fudge_1 = lens_scale_fudge_1; 
	var->lens_scale_fudge_2=lens_scale_fudge_2; 
	var->lens_rad_x = lens_rad_x;
	var->lens_rad_y = lens_rad_y;
	var->lens_rad = lens_rad;
	var->lens_count;	
}
