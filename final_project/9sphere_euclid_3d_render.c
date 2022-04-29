

#include "9sphere_render.c"

int main()
{

  phi = euclid_phi;
  phi_grad = euclid_phi_grad;

  id = create_new_xwd_map(SCREEN_WIDTH, SCREEN_HEIGHT);
  if(id == -1) { printf("failure\n"); exit(0); }

  res = 10;

  do_3d() ;
}
