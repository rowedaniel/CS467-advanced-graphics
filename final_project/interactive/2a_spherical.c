#include "3sphere_morefreedom.c"

int main()
{

  phi = sphere_phi;
  phi_grad = sphere_phi_grad;

  G_init_graphics(SCREEN_WIDTH,SCREEN_HEIGHT);

  delta_t = 0.002;
  max_ray_distance = 300000;

  res = 10;



  do_2d() ;
}
