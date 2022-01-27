


# include "FPToolkit.c"


int main(int argc, char *argv[])
{
	
	// arguments:
	int screen_width, screen_height;
	int framecount;
	double growth;
	
	// other vars:
	double centerx, centery;
	double rad;
	char filename[100];
	double color[3];
	int i;

	if(argc != 5) {
		printf("Usage: pgm_name screen_width screen_height num_frames delta_radius\n");
		return 0;
	}


	screen_width = atoi(argv[1]); screen_height = atoi(argv[2]);
	framecount = atoi(argv[3]);
	growth = atof(argv[4]);

	centerx = screen_width / 2.0;
	centery = screen_height / 2.0;
	rad = 10;
	color[0] = 1;
	color[1] = 0;
	color[2] = 0;

	G_init_graphics(screen_width, screen_height);


	for(i = 0; i < framecount; ++i) {
		G_rgb(0, 0, 0);
		G_clear();

		G_rgb(color[0], color[1], color[2]);
		G_fill_circle(centerx, centery, rad);
		rad += growth;
		color[1] += 1.0 / framecount;

		sprintf(filename, "growingsun%04d.xwd", i);
		G_save_image_to_file(filename);

		G_wait_key();
	}


}
