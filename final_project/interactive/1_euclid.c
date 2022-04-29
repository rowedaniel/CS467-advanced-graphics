
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

  res = 10;
  delta_t = 0.01;
  max_ray_distance = 1000;

  do_2d() ;
}
