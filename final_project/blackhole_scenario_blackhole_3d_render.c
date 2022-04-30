
#include "blackhole_scenario_render.c"

int main(int argc, char * argv[])
{

  phi = blackhole_phi;
  phi_grad = blackhole_phi_grad;

  id = create_new_xwd_map(SCREEN_WIDTH, SCREEN_HEIGHT);
  if(id == -1) { printf("failure\n"); exit(0); }

  delta_t = 0.03;
  max_ray_distance = 10000;

  res = 1;
  y_grid = 0;

  move_per_frame = 0.2;

  int currnum = atoi(argv[1]);
  int countby = atoi(argv[2]);
  int maxnum  = atoi(argv[3]);
  for(fnum = currnum; fnum < maxnum; fnum += countby) {
    do_3d() ;
  }
}
