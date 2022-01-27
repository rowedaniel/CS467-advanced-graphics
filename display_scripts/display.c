#include "FPToolkit.c"

int main(int argc, char *argv[])
{
	G_init_graphics(800, 800);

	G_get_image_from_file(argv[1], 0, 0);
	G_wait_key();

}
