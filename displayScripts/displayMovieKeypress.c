#include "FPToolkit.c"

int main(int argc, char *argv[])
{
	// arguments:
	int screen_width, screen_height;
	int start_frame, end_frame;
	
	// other vars:
	char filename[100];
	int i;

	if(argc != 6) {
		printf("Usage: pgm_name screen_width screen_height prefix_name start_frame end_frame\n");
		return 0;
	}


	screen_width = atoi(argv[1]); screen_height = atoi(argv[2]);
	start_frame  = atoi(argv[4]); end_frame    = atoi(argv[5]);

	G_init_graphics(screen_width, screen_height);


	i = start_frame;
	do {
		G_rgb(0, 0, 0);
		G_clear();
		
		sprintf(filename, "%s%04d.xwd", argv[3], i);
		G_get_image_from_file(filename, 0, 0);

		i += 1;
		if(i > end_frame) {
			i = start_frame;
		}

	} while(G_wait_key() != 'q');


}
