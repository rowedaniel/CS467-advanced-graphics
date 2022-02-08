#include <stdio.h>
#include "FPToolkit.c"
#include "time.h"
#include "M3d_matrix_tools.c"

double sgn(double t) {
	if(t > 0)  { return 1;}
	if(t == 0) { return 0;}
	if(t < 0)  { return -1;}
}


int make_graph(double x[], double y[], int n,
	  double start, double end,
	  double (*X)(double t), double (*Y)(double t))
{
	const double step = 1.0 * (end-start)/n;

	for(int i = 0; i < n; ++i) {
		x[i] = X(start + step*i); 
		y[i] = Y(start + step*i); 
	}
}

int make_graph_step(double x[], double y[],
	 	    double start, double end, double step,
	  	    double (*X)(double t), double (*Y)(double t))
{
	const int n = (end-start)/step + 1;
	make_graph(x,y,n, start,end, X,Y);
	return n;
}

int draw_graph(double x[], double y[], int n) {
	for(int i=0; i<n-1; ++i) {
		G_line(x[i], y[i], x[(i+1)%n], y[(i+1)%n]);
	}
}





int main()
{


	// graph vars
	int n = 10000;
	double x[n], y[n];
	double z[n]; // because the transformation matrix wants an extra 'z'

	// transformation matrix vars
	int Tn = 0;
	int Ttypelist[100] ;
	double Tvlist[100] ;
	double T[4][4], Ti[4][4];

	// set up graphics
	G_init_graphics(800, 800);
	G_rgb(0,0,0);
	G_clear();



	
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// circle :
	//   x^2 + y^2 = 1
	//
	//   u in [0,2*M_PI] gives the entire circle :
	//   x = cos(u) ;  y = sin(u) ; 
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   50.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =  100.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  300.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  500.0 ; Tn++ ;
	// Graph only the part with u in [0.25*M_PI, 1.5*M_PI]
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) { return cos(u); }
		double Y(double u) { return sin(u); }
		make_graph(x, y, n,
			   0.25*M_PI, 1.5*M_PI,
			   &X, &Y);
		G_rgb(1,0,0);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// sum4 :
	//   x^4 + y^4 = 1
	//
	//   u in [-1,1] gives the upper branch of the curve :
	//   x = u  ;  y = pow(1 - u*u*u*u, 0.25) ; 
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   30.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   30.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  170.0 ; Tn++ ;
	// Graph the entire curve.
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) { return u; }
		double Y(double u) { return pow(1 - u*u*u*u, 0.25); }
		make_graph(x, y, n,
			   -1, 1,
			   &X, &Y);
		G_rgb(0,1,0);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// square :
	//   |x| + |y| = 1
	//
	//   u in [0,4] gives the entire square :
	//     u in [0,1],  w = u     ;  x = 1 - w ;  y = w     ; 
	//     u in [1,2],  w = u - 1 ;  x = -w    ;  y = 1 - w ;    
	//     u in [2,3],  w = u - 2 ;  x = w - 1 ;  y = -w    ;
	//     u in [3,4],  w = u - 3 ;  x = w     ;  y = w - 1 ;
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  460.0 ; Tn++ ;
	// Graph the entire curve.
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) {
			if(0.0<=u && u<1.0) { return 1-u; }
			if(1.0<=u && u<2.0) { return 1-u; }
			if(2.0<=u && u<3.0) { return u-3; }
			if(3.0<=u && u<4.0) { return u-3; }
		}
		double Y(double u) {
			if(0<=u && u<1) { return u; }
			if(1<=u && u<2) { return 2-u; }
			if(2<=u && u<3) { return 2-u; }
			if(3<=u && u<4) { return u-4; }
		}
		make_graph(x, y, n,
			   0, 4,
			   &X, &Y);
		G_rgb(0,0,1);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// square (alternate parameterization) :
	//   |x| + |y| = 1
	//
	//   u in [0,2*M_PI] gives the entire square :
	//   x = sgn(cos(u))*[cos(u)]^2 ;   y = sgn(sin(u))*[sin(u)]^2  ;
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  500.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  670.0 ; Tn++ ;
	// Graph the entire curve.
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) { return sgn(cos(u)) * cos(u)*cos(u); }
		double Y(double u) { return sgn(sin(u)) * sin(u)*sin(u); }
		make_graph(x, y, n,
			   0, 2*M_PI,
			   &X, &Y);
		G_rgb(0.4,0.4,0);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// astroid :
	//   sqrt(|x|) + sqrt(|y|) = 1
	//
	//   u in [0,2*M_PI] gives the entire astroid :
	//   x = sgn(cos(u))*[cos(u)]^4 ;  y = sgn(sin(u))*[sin(u)]^4  ;
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   80.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   40.0 ; Tn++ ;
	Ttypelist[Tn] = RZ ; Tvlist[Tn] =   45.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  130.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  650.0 ; Tn++ ;
	// Graph the entire curve.
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) { return sgn(cos(u)) * pow(cos(u), 4); }
		double Y(double u) { return sgn(sin(u)) * pow(sin(u), 4); }
		make_graph(x, y, n,
			   0, 2*M_PI,
			   &X, &Y);
		G_rgb(0.2,0,0.5);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// hyperbola :
	//   x^2 - y^2 = 1 
	//
	//   u in (-inf,inf) gives the right branch :
	//   x = cosh(u) ;   y = sinh(u) ;
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =   70.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   70.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  250.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  150.0 ; Tn++ ;
	// Graph only the part with u in [-1, 1.5]  
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) { return cosh(u); }
		double Y(double u) { return sinh(u); }
		make_graph(x, y, n,
			   -1, 1.5,
			   &X, &Y);
		G_rgb(0,1,0.5);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// parabola :
	//   y = x^2
	//
	//   u in (inf,inf) gives the entire parabola : 
	//   x = u ;  y  = u*u ;
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =  150.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =   50.0 ; Tn++ ;
	Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  140.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  200.0 ; Tn++ ;
	// Graph only the part with u in [-1, 2]
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) { return u; }
		double Y(double u) { return u*u; }
		make_graph(x, y, n,
			   -1, 2,
			   &X, &Y);
		G_rgb(0.5,0.8,0.9);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}
	
	
	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// lemon :
	//   x^2 - (1 - y^2)^3 = 0
	//
	//   u in [0,2*M_PI] gives the entire lemon :
	//   x = [cos(u)]^3  ;   y = sin(u) ;
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =  125.0 ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =  125.0 ; Tn++ ;
	Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  620.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  210.0 ; Tn++ ;
	// Graph the entire curve.
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) { return pow(cos(u), 3); }
		double Y(double u) { return sin(u); }
		make_graph(x, y, n,
			   0, 2*M_PI,
			   &X, &Y);
		G_rgb(0.5,0.2,0.2);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}
	

	//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	// curvething :
	//   x = t-sin(t), y=-cos(t) ;
	Tn = 0 ; 
	Ttypelist[Tn] = SX ; Tvlist[Tn] =  20.0  ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =  20.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  100.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =   30.0 ; Tn++ ;
	// Graph the curve in [0, 6*M_PI].
	{
		M3d_make_movement_sequence_matrix(T, Ti,   Tn,Ttypelist, Tvlist);
		double X(double u) { return u-sin(u); }
		double Y(double u) { return -cos(u); }
		make_graph(x, y, n,
			   0, 6*M_PI,
			   &X, &Y);
		G_rgb(1,1,1);
		M3d_mat_mult_points(x, y, z,
			          T,
			  	  x, y, z, n);	  
		draw_graph(x, y, n);
	}




	G_wait_key();

}
