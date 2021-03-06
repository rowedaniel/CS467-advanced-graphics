

#include "3sphere.c"

int main()
{

  phi = sphere_phi;
  phi_grad = sphere_phi_grad;

  G_init_graphics(SCREEN_WIDTH,SCREEN_HEIGHT);

  light_in_world_space[0] = 2.0;
  light_in_world_space[1] = 3;
  light_in_world_space[2] = 0;

  max_ray_distance = 60000;

  light_resolution_u = 100;
  light_resolution_v = 2;

  do_2d() ;
}
