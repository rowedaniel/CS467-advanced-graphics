#include "blackhole_scenario.c"

int main()
{

  phi = blackhole_phi;
  phi_grad = blackhole_phi_grad;

  G_init_graphics(SCREEN_WIDTH,SCREEN_HEIGHT);

  delta_t = 0.003;
  max_ray_distance = 20000;

  res = 5;
  y_grid = 0;



  do_3d() ;
}
