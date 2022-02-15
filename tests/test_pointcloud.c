#include "graph3D.c"

#define WIN_WIDTH 800
#define WIN_HEIGHT 800

int main()
{

	double zbuff[800][800];
	double renderTransform[4][4];
	double screenTransform[4][4];

	// graph vars
	int n = 10000;
	double x[n], y[n], z[n];

	// transformation matrix vars
	int Tn = 0;
	int Ttypelist[100] ;
	double Tvlist[100] ;
	double T[4][4], Ti[4][4];

	// set up graphics
	G_init_graphics(WIN_WIDTH, WIN_HEIGHT);
	G_rgb(0,0,0);
	G_clear();

	// set up transformation matrix for rendering to screen
	make_win_transform(WIN_WIDTH, WIN_HEIGHT);
	



	
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// cylinder :
	//   x^2 + y^2 = 1
	//
	//   u in [0,2*M_PI], v in [0, 10] :
	//   x = cos(u) ;  y = sin(u) ; z = v
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =  100.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   50.0 ; Tn++ ;
	Ttypelist[Tn] = RY ; Tvlist[Tn] =   90.0 ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =  800.0 ; Tn++ ;
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u, double v) { return cos(u); }
		double Y(double u, double v) { return sin(u); }
		double Z(double u, double v) { return v; }
		G_rgb(1,0,0);
		make_graph_step(0, 2*M_PI, 0.01,
				-10, 10, 0.01,
			        &X, &Y, &Z,
				T);
	}

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// cylinder :
	//   x^2 + y^2 = 1
	//
	//   u in [0,2*M_PI], v in [0, 10] :
	//   x = cos(u) ;  y = sin(u) ; z = v
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   90.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   90.0 ; Tn++ ;
	Ttypelist[Tn] = SZ ; Tvlist[Tn] =   40.0 ; Tn++ ;
	Ttypelist[Tn] = RX ; Tvlist[Tn] =   90.0 ; Tn++ ;
	Ttypelist[Tn] = TZ ; Tvlist[Tn] =  800.0 ; Tn++ ;
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u, double v) { return cos(u); }
		double Y(double u, double v) { return sin(u); }
		double Z(double u, double v) { return v; }
		G_rgb(0,1,0);
		make_graph_step(0, 2*M_PI, 0.01,
				-10, 10, 0.01,
			        &X, &Y, &Z,
				T);
	}

	
	G_wait_key();

}
