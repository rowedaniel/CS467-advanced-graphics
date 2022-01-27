#include "FPToolkit.c"

int main(int argc, char *argv[])
{
	// arguments:
	int screen_width, screen_height;
	int start_frame, end_frame;
	int frame_time;
	int backAndForth;
	
	// other vars:
	char filename[100];
	int i;
	int direction;

	if(argc != 8) {
		printf("Usage: pgm_name screen_width screen_height prefix_name start_frame end_frame frame_time(miliseconds) goBackAndForth(0=no,1=yes)\n");
		return 0;
	}


	screen_width = atoi(argv[1]); screen_height = atoi(argv[2]);
	start_frame  = atoi(argv[4]); end_frame    = atoi(argv[5]);
	frame_time   = atoi(argv[6]);
	backAndForth = atoi(argv[7]);

	G_init_graphics(screen_width, screen_height);


	i = start_frame;
	direction = 1;
	do {
		G_rgb(0, 0, 0);
		G_clear();
		
		// display image
		sprintf(filename, "%s%04d.xwd", argv[3], i);
		G_get_image_from_file(filename, 0, 0);
		G_display_image();

		// increment file counter
		i += direction;
		if(i > end_frame || i <= start_frame) {
			if(backAndForth) {
				direction *= -1;
			} else {
				i = start_frame;
			}
		}

		// check if user has signalled to stop
		if(G_no_wait_key() == 'q') {
			break;
		}

		// sleep until next frame
		usleep(frame_time);

	} while(1);


}
