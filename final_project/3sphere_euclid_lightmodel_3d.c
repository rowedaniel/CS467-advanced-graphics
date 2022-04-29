#include "3sphere.c"

int main()
{

  phi = euclid_phi;
  phi_grad = euclid_phi_grad;
  /*
  */

  /*
  phi = sphere_phi;
  phi_grad = sphere_phi_grad;
  */

  G_init_graphics(SCREEN_WIDTH,SCREEN_HEIGHT);

  do_lightmodel = 1;

  light_in_world_space[0] = 2.0;
  light_in_world_space[1] = 3;
  light_in_world_space[2] = 0;

  light_resolution_u = 100;
  light_resolution_v = 100;

  do_3d() ;
}
