
#include "blackhole_scenario_render.c"

int main()
{

  phi = blackhole_phi;
  phi_grad = blackhole_phi_grad;

  id = create_new_xwd_map(SCREEN_WIDTH, SCREEN_HEIGHT);
  if(id == -1) { printf("failure\n"); exit(0); }

  delta_t = 0.003;
  max_ray_distance = 20000;

  res = 1;
  y_grid = 0;



  do_3d() ;
}
