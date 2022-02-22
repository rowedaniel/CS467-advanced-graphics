#include <stdio.h>
#include <stdlib.h>
#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "xwd_tools_03.c"


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



int get_normal(double u, double v,
		double ustep, double vstep,
	  	double (*X)(double u, double v),
	  	double (*Y)(double u, double v),
	  	double (*Z)(double u, double v),
		double T[4][4],
		double N[3])
{
	double P[3]  = {X(u,v), Y(u,v), Z(u,v)};
	double P2[3] = {X(u+ustep,v), Y(u+ustep,v), Z(u+ustep,v)};
	double P3[3] = {X(u,v+vstep), Y(u,v+vstep), Z(u,v+vstep)};
	// get eye, P, and N
	// transform points to camera-space
	M3d_mat_mult_pt(P,  T, P);
	M3d_mat_mult_pt(P2, T, P2);
	M3d_mat_mult_pt(P3, T, P3);

	// subtract P2-P, P3-P to get vectors, not points
	M3d_vector_mult_const(P2, P2, -1);
	M3d_vector_add(P2, P2, P);
	M3d_vector_mult_const(P2, P2, -1);
	
	M3d_vector_mult_const(P3, P3, -1);
	M3d_vector_add(P3, P3, P);
	M3d_vector_mult_const(P3, P3, -1);


	M3d_x_product(N, P2, P3);

}


double calc_area_def(double u, double v,
		   double ustep, double vstep,
	  	   double (*X)(double u, double v),
	  	   double (*Y)(double u, double v),
	  	   double (*Z)(double u, double v),
		   double T[4][4]
	       )
{
	double N[3];
	get_normal(u,v,ustep,vstep,X,Y,Z,T,N);
	double mag1 = M3d_magnitude(N);
	return mag1;

}








// not intended to be used outside of this file
int make_graph_base(
	       double ustart, double uend, double ustep,
	       double vstart, double vend, double vstep,
	       double (*X)(double u, double v),
	       double (*Y)(double u, double v),
	       double (*Z)(double u, double v),
	       double T[4][4],
	       int (*light_model_wrapper)(double u, double v,
		       			 double eye[3], double P[3], double N[3],
					 double rgb_out[3])
	       )
{
	double u, v;
	double bigu, bigv;
	double comp_area_def;


	// keep trying random comparison points until one of them isn't 0.
	srand(4242);
	for(int i=0; i<1000; ++i){
		const double u_random = (uend-ustart) * (1.0*rand()/RAND_MAX) + ustart;
		const double v_random = (vend-vstart) * (1.0*rand()/RAND_MAX) + vstart;
		comp_area_def = calc_area_def(u_random, v_random, ustep,vstep, X,Y,Z, T);
		if (comp_area_def == 0.0) {
			break;
		}
	} 
	if (comp_area_def == 0.0) {
		comp_area_def = 1;
	}

	for(bigu = ustart; bigu < uend; bigu += ustep) {


		for(bigv = vstart; bigv < vend; bigv += vstep) {

			// refine to make sure the resolution is the same
			double area_def = calc_area_def(bigu, bigv, ustep,vstep, X,Y,Z, T);
			if(area_def * area_def == area_def) { area_def = 0; } // nan, inf, or 0
			double ratio = area_def/comp_area_def;
			//double ratio = 1;
			if(ratio != 1) {
				//printf("ratio: %i, area_def:%lf, comp_area_def:%lf\n", ratio, area_def, comp_area_def);
			}
			for(double u = bigu; u < bigu+ustep; u += ustep / ratio) {
				for(double v = bigv; v < bigv+vstep; v += vstep / ratio) {

	
	
					double P[3] = {X(u,v), Y(u,v), Z(u,v)};
					
					// first, apply the given transformation to move into camera-space
					M3d_mat_mult_pt(P, T, P);
		
					// don't bother rendering if point is behind camera
					if(P[2] <= 0) { continue; } 		
					// also, don't bother rendering if point is not in cone of vision
					if(fabs(P[0]/P[2]) > 1) { continue; }
					if(fabs(P[1]/P[2]) > 1) { continue; }
		
		
		
					// next, do the light model
					double new_rgb[3];
					double eye[3] = {0, 0, 0};
					double N[3];
		
					get_normal(u, v, ustep, vstep,
						   X, Y, Z,
						   T,
						   N);
		
					light_model_wrapper(u, v,   eye, P, N, new_rgb);
					//Light_Model(inherent_rgb, eye, P, N, new_rgb);
		
		
		
					G_rgb(new_rgb[0], new_rgb[1], new_rgb[2]);
					//G_rgb(inherent_rgb[0], inherent_rgb[1], inherent_rgb[2]);
		
		
					plot(P); 
					//printf("u,v: (%lf, %lf)\n", u, v);
	
				}
			}
	
		}

	}

}



// using inherent rgbs per object:





int make_graph(
	       double ustart, double uend, double ustep,
	       double vstart, double vend, double vstep,
	       double (*X)(double u, double v),
	       double (*Y)(double u, double v),
	       double (*Z)(double u, double v),
	       double T[4][4],
	       double inherent_rgb[3])
{

	int rgb_light_model_wrapper(double u, double v,
		       		   double eye[3], double P[3], double N[3],
				   double new_rgb[3])
	{
		Light_Model(inherent_rgb, eye, P, N, new_rgb);
	}


	make_graph_base(ustart,uend,ustep,
			vstart,vend,vstep,
			X,Y,Z,
			T,
			rgb_light_model_wrapper);
}



int make_graph_image(
	       double ustart, double uend, double ustep,
	       double vstart, double vend, double vstep,
	       double (*X)(double u, double v),
	       double (*Y)(double u, double v),
	       double (*Z)(double u, double v),
	       double T[4][4],
	       int image_ID, int image_width, int image_height)
{

	int rgb_light_model_wrapper(double u, double v,
		       		   double eye[3], double P[3], double N[3],
				   double new_rgb[3])
	{
		double rgb[3];
		int x = floor( image_width * (u-ustart) / (uend-ustart) );
		int y = floor( image_height* (v-vstart) / (vend-vstart) );
		int e = get_xwd_map_color(image_ID,  x,y,rgb);
		if(e == -1) { return -1; }
		Light_Model(rgb, eye, P, N, new_rgb);
	}


	make_graph_base(ustart,uend,ustep,
			vstart,vend,vstep,
			X,Y,Z,
			T,
			rgb_light_model_wrapper);
}


