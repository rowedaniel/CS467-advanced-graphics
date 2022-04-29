

#include "9sphere_render.c"

int main()
{

  phi = sphere_phi;
  phi_grad = sphere_phi_grad;

  id = create_new_xwd_map(SCREEN_WIDTH, SCREEN_HEIGHT);
  if(id == -1) { printf("failure\n"); exit(0); }

  do_lightmodel = 1;
  res = 5;

  light_in_world_space[0] = 2.0;
  light_in_world_space[1] = 3;
  light_in_world_space[2] = 0;

  light_resolution_u = 100;
  light_resolution_v = 100;

  do_3d() ;
}
