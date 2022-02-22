#include "graph3D.c"



// declare global vars

#define M 10
double numobjects = 0;
double inherent_rgb[M][3] ;







void draw_all_objects(double V[4][4], double T[M][4][4],
		double (*X[M]) (double u, double v),
		double (*Y[M]) (double u, double v),
		double (*Z[M]) (double u, double v),
		double uStart[M], double uEnd[M], double uStep[M],
		double vStart[M], double vEnd[M], double vStep[M]
		)
{
  double t[4][4];

  for(int onum = 0; onum < numobjects; ++onum) {

    M3d_mat_mult(t, V, T[onum]);
    make_graph(uStart[onum], uEnd[onum], uStep[onum],
		    vStart[onum], vEnd[onum], vStep[onum],
    	            (X[onum]), (Y[onum]), (Z[onum]), 
    		    t,
		    inherent_rgb[onum]);
  }
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


int main ()
{
  int onum;

  // stuff for moving each object to the right location
  double V[M][4][4], Vi[M][4][4] ;
  int nl ;
  int tlist[100] ;
  double plist[100] ;
  // array of generation functions for each shape
  double (*X[M]) (double u, double v);
  double (*Y[M]) (double u, double v);
  double (*Z[M]) (double u, double v);

  // array of start/stop/step values for each parameterization
  double uStart[M], uEnd[M], uStep[M];
  double vStart[M], vEnd[M], vStep[M];



  // initialize pointCloud stuff
  make_win_transform(600, 600);
  // init graphics
  G_init_graphics (win_width, win_height) ;




  //---------------------------------------------------------------------
  // Design the model by moving objects from their own OBJECT SPACE 
  // to a common WORLD SPACE :
  //---------------------------------------------------------------------

  onum = 0 ;  // current object number

  // The sphere has a radius of 1.0, centered at the origin.
  double sphereX(double u, double v) { return sqrt(1-v*v)*cos(u); }
  double sphereY(double u, double v) { return v; }
  double sphereZ(double u, double v) { return sqrt(1-v*v)*sin(u); }
  
  // The cylinder has radius of 0.1, centered at origin, 
  // lying on the x-axis from x = -1 to x = 1.
  double cylinderX(double u, double v) { return u; }
  double cylinderY(double u, double v) { return 0.1*cos(v); }
  double cylinderZ(double u, double v) { return 0.1*sin(v); }




  // Build an origin point by placing a sphere at 0,0,0 
  inherent_rgb[onum][0] = 1.0 ;
  inherent_rgb[onum][1] = 0.8 ;
  inherent_rgb[onum][2] = 0.0 ;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.25 ; nl++ ;  
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;

  X[onum] = sphereX; 
  Y[onum] = sphereY; 
  Z[onum] = sphereZ; 

  uStart[onum] = 0; 		uEnd[onum] = 2*M_PI,	 uStep[onum] = 0.01;
  vStart[onum] = -1;		vEnd[onum] = 1,		 vStep[onum] = 0.01;

  onum++ ;


  // Build a +x axis with a cylinder.
  inherent_rgb[onum][0] = 1.0 ;
  inherent_rgb[onum][1] = 0.2 ;
  inherent_rgb[onum][2] = 0.2 ;

  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 2.00 ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  X[onum] = cylinderX; 
  Y[onum] = cylinderY; 
  Z[onum] = cylinderZ; 

  uStart[onum] = -1; 		uEnd[onum] = 1,	 	 uStep[onum] = 0.001;
  vStart[onum] =  0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.001;

  onum++ ;


  // Build a +y axis with the cylinder file.
  inherent_rgb[onum][0] = 0.3 ;
  inherent_rgb[onum][1] = 1.0 ;
  inherent_rgb[onum][2] = 0.2 ;
  
  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 2.00 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] =  90  ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  X[onum] = cylinderX; 
  Y[onum] = cylinderY; 
  Z[onum] = cylinderZ; 

  uStart[onum] = -1; 		uEnd[onum] = 1,	 	 uStep[onum] = 0.001;
  vStart[onum] =  0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.001;

  onum++ ;


  // Build a +z axis with a cylinder.
  inherent_rgb[onum][0] = 0.3 ;
  inherent_rgb[onum][1] = 0.2 ;
  inherent_rgb[onum][2] = 1.0 ;

  nl = 0 ;
  tlist[nl] = TX ; plist[nl] = 1.00 ; nl++ ;
  tlist[nl] = SX ; plist[nl] = 2.00 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -90  ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  X[onum] = cylinderX; 
  Y[onum] = cylinderY; 
  Z[onum] = cylinderZ; 

  uStart[onum] = -1; 		uEnd[onum] = 1,	 	 uStep[onum] = 0.001;
  vStart[onum] =  0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.001;

  onum++ ;


  // Build a diagonal cylinder going from +x to +z 
  inherent_rgb[onum][0] = 0.5 ;
  inherent_rgb[onum][1] = 1.0 ;
  inherent_rgb[onum][2] = 0.4 ;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 2.0*sqrt(2.0); nl++ ;
  tlist[nl] = TZ ; plist[nl] = 3.00 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 45   ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  X[onum] = cylinderX; 
  Y[onum] = cylinderY; 
  Z[onum] = cylinderZ; 

  uStart[onum] = -1; 		uEnd[onum] = 1,	 	 uStep[onum] = 0.001;
  vStart[onum] =  0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.001;

  onum++ ;


  // Build a diagonal cylinder going from +x to +z 
  inherent_rgb[onum][0] = 0.5 ;
  inherent_rgb[onum][1] = 1.0 ;
  inherent_rgb[onum][2] = 0.4 ;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 2.0*sqrt(2.0); nl++ ;
  tlist[nl] = TY ; plist[nl] = 3.00 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = -45  ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  X[onum] = cylinderX; 
  Y[onum] = cylinderY; 
  Z[onum] = cylinderZ; 

  uStart[onum] = -1; 		uEnd[onum] = 1,	 	 uStep[onum] = 0.001;
  vStart[onum] =  0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.001;

  onum++ ;


  // Build a diagonal cylinder going from +y to +z
  inherent_rgb[onum][0] = 0.5 ;
  inherent_rgb[onum][1] = 1.0 ;
  inherent_rgb[onum][2] = 0.4 ;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 2.0*sqrt(2.0); nl++ ;
  tlist[nl] = RY ; plist[nl] = 90   ; nl++ ;
  tlist[nl] = TY ; plist[nl] = 3.00 ; nl++ ;
  tlist[nl] = RX ; plist[nl] = 45   ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  X[onum] = cylinderX; 
  Y[onum] = cylinderY; 
  Z[onum] = cylinderZ; 

  uStart[onum] = -1; 		uEnd[onum] = 1,	 	 uStep[onum] = 0.001;
  vStart[onum] =  0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.001;

  onum++ ;


  // Build a corner point for +x out of the sphere
  inherent_rgb[onum][0] = 0.8 ;
  inherent_rgb[onum][1] = 1.0 ;
  inherent_rgb[onum][2] = 0.7 ;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.25 ; nl++ ;  
  tlist[nl] = TX ; plist[nl] = 4.20 ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;

  X[onum] = sphereX; 
  Y[onum] = sphereY; 
  Z[onum] = sphereZ; 

  uStart[onum] = 0; 		uEnd[onum] = 2*M_PI,	 uStep[onum] = 0.001;
  vStart[onum] = -1;		vEnd[onum] = 1,		 vStep[onum] = 0.001;

  onum++ ;


  // Build an corner point for +z out of the sphere file.
  inherent_rgb[onum][0] = 0.8 ;
  inherent_rgb[onum][1] = 1.0 ;
  inherent_rgb[onum][2] = 0.7 ;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.25 ; nl++ ;  
  tlist[nl] = TZ ; plist[nl] = 4.20 ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;

  X[onum] = sphereX; 
  Y[onum] = sphereY; 
  Z[onum] = sphereZ; 

  uStart[onum] = 0; 		uEnd[onum] = 2*M_PI,	 uStep[onum] = 0.001;
  vStart[onum] = -1;		vEnd[onum] = 1,		 vStep[onum] = 0.001;

  onum++ ;


  // Build an corner point for +y out of the sphere file.
  inherent_rgb[onum][0] = 0.8 ;
  inherent_rgb[onum][1] = 1.0 ;
  inherent_rgb[onum][2] = 0.7 ;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 0.25 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 0.25 ; nl++ ;  
  tlist[nl] = TY ; plist[nl] = 4.20 ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;

  X[onum] = sphereX; 
  Y[onum] = sphereY; 
  Z[onum] = sphereZ; 

  uStart[onum] = 0; 		uEnd[onum] = 2*M_PI,	 uStep[onum] = 0.001;
  vStart[onum] = -1;		vEnd[onum] = 1,		 vStep[onum] = 0.001;

  onum++ ;






  //---------------------------------------------------------------------
  numobjects = onum ;
  //---------------------------------------------------------------------


  // eye stuff
  double eye[3], coi[3], up[3] ;
  double v[4][4], vi[4][4];
  int i;

  int fnum ;
  int max_fnum;
  double t ;

  fnum = 0 ;
  max_fnum = 300;


  // light model setup
  double light_in_world_space[3] = {1, 2, 0};
  
  AMBIENT = 0.2 ;
  MAX_DIFFUSE = 0.5 ;
  SPECPOW = 30 ;

  // file saving setup
  char filename[100];



  while (1) {

    t = 0.01*fnum ;

    eye[0] = 15*cos(2*M_PI*t) ; 
    eye[1] =  6*t ; 
    eye[2] =  7*sin(2*M_PI*t) ; 

    // printf("t = %lf   eye = %lf %lf %lf\n",t, eye[0],eye[1],eye[2]) ;

    coi[0] =  1.0 ;
    coi[1] =  2.0 ; 
    coi[2] =  0.5 ;

    up[0]  = eye[0] ; 
    up[1]  = eye[1] + 1 ;
    up[2]  = eye[2] ; 

    //------------------------------- put your code here!!!!!!!!!!!!

    // view the screen
    G_rgb(0,0,0);
    G_clear();
    printf("rendering frame %d\n", fnum);


    // render
    M3d_view(v,vi, eye, coi, up);

    M3d_mat_mult_pt(light_in_eye_space, v, light_in_world_space);
    
    make_z_buff();
    draw_all_objects(v, V,
		     X, Y, Z,
		     uStart, uEnd, uStep,
		     vStart, vEnd, vStep
		     );






    G_display_image();

    sprintf(filename, "growingsun%04d.xwd", fnum);
    G_save_image_to_file(filename);


    if(G_no_wait_key() == 'q' || fnum >= max_fnum) {
      break;
    }


    fnum++ ;
  } // end while (1)


  while ('q' != G_wait_key()) {}


}

