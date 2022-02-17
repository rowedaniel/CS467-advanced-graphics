#include <stdio.h>
#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


#define MAX_WIN_WIDTH 2000
#define MAX_WIN_HEIGHT 2000


double win_width, win_height;
double zbuff[MAX_WIN_WIDTH][MAX_WIN_HEIGHT];
double screenTransform[4][4];

int make_win_transform(double w, double h) {
	
	int Tn = 0;
	int Ttypelist[100] ;
	double Tvlist[100] ;
	double Ti[4][4];

	Ttypelist[Tn] = SX ; Tvlist[Tn] =  w / 2.0  ; Tn++ ;
	Ttypelist[Tn] = SY ; Tvlist[Tn] =  h / 2.0 ; Tn++ ;
	Ttypelist[Tn] = TX ; Tvlist[Tn] =  w / 2.0 ; Tn++ ;
	Ttypelist[Tn] = TY ; Tvlist[Tn] =  h / 2.0 ; Tn++ ;
	M3d_make_movement_sequence_matrix(screenTransform, Ti,   Tn,Ttypelist, Tvlist);

	win_width = w;
	win_height = h;
}

int make_z_buff() {
	for(int i=0; i< win_width; ++i) {
		for(int j=0; j<win_height; ++j) {
			zbuff[i][j] = 9999999;
		}
	}
}


int plot(double x, double y, double z, double T[4][4])
{
	double P[3] = {x, y, z};
	
	// first, apply the given transformation to move into camera-space
	M3d_mat_mult_pt(P, T, P);
	if(P[2] <= 0) { return 0; } // don't bother rendering if point is behind camera

	// scale to z=1 plane
	P[0] /= P[2];
	P[1] /= P[2];
	double zbar = P[2];

	// move to screen-space
	M3d_mat_mult_pt(P, screenTransform, P);

	int xbarbar = (int)P[0];
	int ybarbar = (int)P[1];

	// if it's closer to the camera than the previous point, go ahead and render it
	// also make sure it's actually inside the window
	if(zbar < zbuff[xbarbar][ybarbar] && 1) {
	   //0 <= xbarbar && xbarbar < win_width &&
	   //0 <= ybarbar && ybarbar < win_height) {
		G_point(xbarbar, ybarbar);
		zbuff[xbarbar][ybarbar] = (int)P[2];
	}

	return 1;
}

int make_graph_step(double ustart, double uend, double ustep,
	  	    double vstart, double vend, double vstep,
	  	    double (*X)(double u, double v),
	  	    double (*Y)(double u, double v),
	  	    double (*Z)(double u, double v),
		    double T[4][4])
{
	double u, v;
	for(u = ustart; u < uend; u += ustep) {
		for(v = vstart; v < vend; v += vstep) {
			plot(X(u,v), Y(u,v), Z(u,v), T); 
		}
	}

}



