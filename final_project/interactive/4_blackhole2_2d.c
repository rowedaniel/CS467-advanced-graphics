#include "blackhole_scenario.c"

int main()
{

  phi = blackhole_phi;
  phi_grad = blackhole_phi_grad;

  G_init_graphics(SCREEN_WIDTH,SCREEN_HEIGHT);

  delta_t = 0.03;
  max_ray_distance = 10000;

  res = 1;



  do_2d() ;
}
