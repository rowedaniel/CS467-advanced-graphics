


# include "FPToolkit.c"


int main(int argc, char *argv[])
{
	printf("argc: %d\n", argc);
	
	// arguments:
	int swidth, sheight;
	int framecount;
	double growth;
	
	// other vars:
	double centerx, centery;
	double rad;
	char filename[100];
	double color[3];
	int i;

	if(argc != 5) {
		printf("Usage: pgm <width> <height> <frames> <growth>\n");
		return 0;
	}


	swidth = atoi(argv[1]); sheight = atoi(argv[2]);
	framecount = atoi(argv[3]);
	growth = atof(argv[4]);

	centerx = swidth / 2.0;
	centery = sheight / 2.0;
	rad = 10;
	color[0] = 1;
	color[1] = 0;
	color[2] = 0;

	G_init_graphics(swidth, sheight);


	for(i = 0; i < framecount; ++i) {
		G_rgb(0, 0, 0);
		G_clear();

		G_rgb(color[0], color[1], color[2]);
		G_fill_circle(centerx, centery, rad);
		rad += growth;
		color[1] += 1.0 / framecount;

		sprintf(filename, "growingsun%i.xwd", i);
		G_save_image_to_file(filename);

		G_wait_key();
	}


}
