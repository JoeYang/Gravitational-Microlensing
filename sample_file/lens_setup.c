#include<stdio.h>
#include<stdlib.h>
#include<math.h>

int main() {

  /* Two parameters for the lensing system: convergence (kappa) and shear (gamma).
     Convergence can be broken into a smooth component (kappa_c) and a clumpy stellar
     component (kappa_star), such that kappa = kappa_star + kappa_c */

  float kappa_star = 0.394;           // convergence in stars (sometimes written as sigma instead of kappa)
  float kappa_c = 0.;                 // convergence in smooth matter
  float gamma = 0.395;                // shear
  float kappa = kappa_star + kappa_c; // total convergence

  /* Use a single mass for all the microlenses at the moment. This is measured in solar masses. */
  float lens_mass = 1.;    // microlens mass

  /* Source plane (where the light rays are collected) is a square with side length source_scale,
     measured in Einstein Radii. */

  float source_scale = 24.;       // size of the source plane in Einstein Radii

  /* Size of the image plane (from which the light rays are shot) depends on the desired size of the 
     source plane and the lensing parameters. Note that the image plane is rectangular -- 2d ray 
     positions should be randomly generated in the ranges [-image_scale_x, image_scale_x] and
     [-image_scale_y, image_scale_y].
     Fudge factor values were determined empirically -- leaving them set to the following values 
     will allow direct comparison of output with single core tree code */

  float image_scale_fudge = 0.1;  // image plane scale fudge factor
  float image_scale_x = (0.5 + 2.0*image_scale_fudge) * source_scale / (1.0 - kappa - gamma);
  float image_scale_y = (0.5 + 2.0*image_scale_fudge) * source_scale / (1.0 - kappa + gamma);

  /* Radius of the lens plane (where the microlenses are distributed) also depends on the desired size
     of the source plane and the lensing parameters. We use a circular lens plane of radius lens_rad.
     Some notes apply about the fudge factors. */

  float lens_scale_fudge_1 = 5.0; // lens plane scale fudge factor 1
  float lens_scale_fudge_2 = 2.0; // lens plane scale fudge factor 2
  float lens_rad_x = (0.5*source_scale + lens_scale_fudge_1) / (1.0 - kappa + gamma);
  float lens_rad_y = (0.5*source_scale + lens_scale_fudge_1) / (1.0 - kappa - gamma);
  float lens_rad = sqrt(lens_rad_x*lens_rad_x + lens_rad_y*lens_rad_y) + lens_scale_fudge_2;

  /* lens_count is the required number of microlenses. These should be uniformly randomly distributed 
     within a circle of radius lens_rad. */  

  size_t lens_count = (size_t)(kappa_star * lens_rad*lens_rad / lens_mass + 0.5); // TODO: Eventually this should be the avg. lens mass
  
}
