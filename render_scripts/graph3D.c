#include <stdio.h>
#include "FPToolkit.c"
#include "M3d_matrix_tools.c"


#define MAX_WIN_WIDTH 2000
#define MAX_WIN_HEIGHT 2000


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// To support the light model :
double light_in_eye_space[3] ;
double AMBIENT      = 0.2 ;
double MAX_DIFFUSE  = 0.5 ;
double SPECPOW      = 50 ;



int Light_Model (double irgb[3],
                 double s[3],
                 double p[3],
                 double n[3],
                 double argb[3])
// s,p,n in eyespace

// irgb == inherent color of object (input to this function)
// s = location of start of ray (probably the eye)
// p = point on object (input to this function)
// n = normal to the object at p (input to this function)
// argb == actual color of object (output of this function)
// globals : AMBIENT, MAX_DIFFUSE, SPECPOW, light_in_eye_space[3]

// return 1 if successful, 0 if error
{

  double len ;
  double N[3] ; 
  len = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]) ;
  if (len == 0) return 0 ;
  N[0] = n[0]/len ;  N[1] = n[1]/len ;  N[2] = n[2]/len ;

  double E[3] ;
  E[0] = s[0] - p[0] ; 
  E[1] = s[1] - p[1] ; 
  E[2] = s[2] - p[2] ; 
  len = sqrt(E[0]*E[0] + E[1]*E[1] + E[2]*E[2]) ;
  if (len == 0) return 0 ;
  E[0] /= len ;  E[1] /= len ;  E[2] /= len ;
  double NdotE = N[0]*E[0] + N[1]*E[1] + N[2]*E[2] ;

  double L[3] ;
  L[0] = light_in_eye_space[0] - p[0] ; 
  L[1] = light_in_eye_space[1] - p[1] ; 
  L[2] = light_in_eye_space[2] - p[2] ; 
  len = sqrt(L[0]*L[0] + L[1]*L[1] + L[2]*L[2]) ;
  if (len == 0) return 0 ;
  L[0] /= len ;  L[1] /= len ;  L[2] /= len ;
  double NdotL = N[0]*L[0] + N[1]*L[1] + N[2]*L[2] ;





  double max_ambient_and_diffuse = AMBIENT + MAX_DIFFUSE ;
     // this needs to occur BEFORE you possibly jump to LLL below




  double intensity ;
  if (NdotL*NdotE < 0) {
    // eye and light are on opposite sides of polygon
    intensity = AMBIENT ; 
    goto LLL ;
  } else if ((NdotL < 0) && (NdotE < 0)) {
    // eye and light on same side but normal pointing "wrong" way
    N[0] *= (-1.0) ;    N[1] *= (-1.0) ;    N[2] *= (-1.0) ; 
    NdotL *= (-1.0) ;
    NdotE *= (-1.0) ;   // don't use NdotE below, probably should eliminate this
  }


  // ignore Blinn's variant
  double R[3] ; // Reflection vector of incoming light
  R[0] = 2*NdotL*N[0] - L[0] ;
  R[1] = 2*NdotL*N[1] - L[1] ;
  R[2] = 2*NdotL*N[2] - L[2] ;

  double EdotR = E[0]*R[0] + E[1]*R[1] + E[2]*R[2] ;

  double diffuse ;
  if (NdotL <= 0.0) { diffuse = 0.0 ; }
  else { diffuse = MAX_DIFFUSE*NdotL ; }

  double specular ;
  if (EdotR <= 0.0) { specular = 0.0 ; }
  else { specular = (1.0 - max_ambient_and_diffuse)*pow(EdotR,SPECPOW) ;}

  // printf("%lf %lf\n",diffuse,specular) ;
  intensity = AMBIENT + diffuse + specular ;



 LLL : ;

  double f,g ;
  if (intensity <= max_ambient_and_diffuse) {
    f = intensity / max_ambient_and_diffuse ;
    argb[0] = f * irgb[0] ;
    argb[1] = f * irgb[1] ;
    argb[2] = f * irgb[2] ;
  } else {
    f = (intensity - max_ambient_and_diffuse) / 
                           (1.0 - max_ambient_and_diffuse) ;
    g = 1.0 - f ;
    argb[0] = g * irgb[0] + f ;
    argb[1] = g * irgb[1] + f ;
    argb[2] = g * irgb[2] + f ;
  }

  return 1 ;
}






void light_model (double irgb[3],
                  double *P, double *P2, double *P3, 
                  double argb[3])
// irgb == inherent color of object (input to this function)
// xx[],yy[],zz[] are points in the polygon
// argb == actual color of object (output of this function)
{
  double Eye[3] ;
  Eye[0] = 0 ; Eye[1] = 0 ; Eye[2] = 0 ; 

  double a[3] ;
  M3d_vector_mult_const(a, P, -1);
  M3d_vector_add(a, a, P2);

  double b[3] ;
  M3d_vector_mult_const(b, P, -1);
  M3d_vector_add(b, b, P3);
 
  double N[3] ;
  M3d_x_product (N, a,b) ;

  Light_Model (irgb, Eye, P, N, argb) ;
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++






// to support graphing
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


int plot(double P[3])
{



	


	// scale to z=1 plane
	P[0] /= P[2];
	P[1] /= P[2];
	double zbar = P[2];

	// move to screen-space
	M3d_mat_mult_pt(P, screenTransform, P);

	int xbarbar = floor(P[0]);
	int ybarbar = floor(P[1]);

	// if it's closer to the camera than the previous point, go ahead and render it
	// also make sure it's actually inside the window
	if(zbar < zbuff[xbarbar][ybarbar] && 1) {
	   //0 <= xbarbar && xbarbar < win_width &&
	   //0 <= ybarbar && ybarbar < win_height) {
	   	

		G_point(xbarbar, ybarbar);
		zbuff[xbarbar][ybarbar] = zbar;
	}

	return 1;
}

int make_graph_step(double ustart, double uend, double ustep,
	  	    double vstart, double vend, double vstep,
	  	    double (*X)(double u, double v),
	  	    double (*Y)(double u, double v),
	  	    double (*Z)(double u, double v),
		    double T[4][4],
		    double inherent_rgb[3])
{
	double u, v;
	for(u = ustart; u < uend; u += ustep) {
		for(v = vstart; v < vend; v += vstep) {
			double P[3] = {X(u,v), Y(u,v), Z(u,v)};
			
			// first, apply the given transformation to move into camera-space
			M3d_mat_mult_pt(P, T, P);

			// don't bother rendering if point is behind camera
			if(P[2] <= 0) { continue; } 		

			// next, do the lightmodel
			double new_rgb[3];
			double P2[3] = {X(u+ustep,v), Y(u+ustep,v), Z(u+ustep,v)};
			double P3[3] = {X(u,v+vstep), Y(u,v+vstep), Z(u,v+vstep)};
			// transform points to camera-space
			M3d_mat_mult_pt(P2, T, P2);
			M3d_mat_mult_pt(P3, T, P3);
			// do light model
			light_model(inherent_rgb, P, P2, P3, new_rgb);
			G_rgb(new_rgb[0], new_rgb[1], new_rgb[2]);
			//G_rgb(inherent_rgb[0], inherent_rgb[1], inherent_rgb[2]);


			plot(P); 
		}
	}

}



