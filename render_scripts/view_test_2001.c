#include "graph3D.c"



// declare global vars

#define M 100
double numobjects = 0;
int image_IDs[M];
int image_size[M][2];
double inherent_rgb[M][3];







void draw_all_objects(double V[4][4], double T[M][4][4],
		int (*f[M]) (double u, double v, double P[3]),
		double uStart[M], double uEnd[M], double uStep[M],
		double vStart[M], double vEnd[M], double vStep[M]
		)
{
  double t[4][4];

  for(int onum = 0; onum < numobjects; ++onum) {

    M3d_mat_mult(t, V, T[onum]);
    make_graph_inherentrgb_1func(uStart[onum], uEnd[onum], uStep[onum],
	       vStart[onum], vEnd[onum], vStep[onum],
    	       f[onum], 
               t,
	       inherent_rgb[onum]
	       );

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
  int (*f[M]) (double u, double v, double P[3]);

  // array of start/stop/step values for each parameterization
  double uStart[M], uEnd[M], uStep[M];
  double vStart[M], vEnd[M], vStep[M];

  // stuff for loading image files
  int error;



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
  double sphere(double u, double v, double P[3]) {
	  const double a = 1-v*v;
	  if(a<0) { return 0; }
	  const double b = sqrt(a);
  	  P[0] = b*cos(u);
      	  P[1] =  v;
      	  P[2] = b*sin(u);
    	  return 1;
  }

  // Torus has center radius 1, outer radius 0.1
  int torus(double u, double v, double P[3]) {
	  const double torus_centerR = 1;
	  const double torus_outerR = 0.1;
	  const double torus_squareness = 1.0 / 4.0;
	  double f(double x) {
	  	double a = cos(x);
	  	return ((a>0)-(a<0)) * pow(fabs(a), torus_squareness);
	  }
	  P[0] = (torus_centerR + torus_outerR*f(v)) * cos(u);
	  P[1] = (torus_centerR + torus_outerR*f(v)) * sin(u);
	  P[2] = torus_outerR*sin(v);
	  return 1;
  }
  
  // The cylinder has radius of 0.1, centered at origin, 
  // lying on the x-axis from x = -1 to x = 1.
  int cylinder(double u, double v, double P[3]) {
	  double rad = 0.1;
  	  P[0] = u;
  	  P[1] =  rad*cos(v);
  	  P[2] =  rad*sin(v);
	  return 1;
  }

  // door bays are circles (filled in), but with a rectangular slot cut out
  int doorbay(double u, double v, double P[3]) {
	  const double doorbay_angle = M_PI/9;
	  const double doorbay_inset = 0.8;
	  const double doorbay_x = doorbay_inset * cos(doorbay_angle);
	  const double doorbay_y = doorbay_inset * sin(doorbay_angle);

	  const double x=v*cos(u);
	  const double y=v*sin(u);

	  if(-doorbay_x < x && x < doorbay_x && -doorbay_y < y && y < doorbay_y) {
	    return 0;
	  }

	  P[0] =  x;
	  P[1] =  y;
	  P[2] =  0;
	  return 1;
  }


  double station_color[3] = {0.8, 0.8, 0.9};
  double station_highlight[3] = {0.7, 0.3, 0.5};



  //image_IDs[onum] = init_xwd_map_from_file("clock.xwd");
  //if(image_IDs[onum] == -1) { printf("File load failure!\n"); exit(0); }
  //error = get_xwd_map_dimensions(image_IDs[onum],image_size[onum]); 
  //if(error == -1) { printf("File load failure!\n"); exit(0); }


  // Build central axis cylinder.
  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 2.00 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 6.00 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 6.00 ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  f[onum] = cylinder; 

  uStart[onum] = -1;		uEnd[onum] = 1, 	 uStep[onum] = 0.01;
  vStart[onum] = 0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.01;

  onum++ ;


  for(int i=0; i<2; ++i) {
	  // Build radial torusi
	  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
	
	  nl = 0 ;
	  tlist[nl] = SX ; plist[nl] = 4.00 ; nl++ ;
	  tlist[nl] = SY ; plist[nl] = 4.00 ; nl++ ;
	  tlist[nl] = SZ ; plist[nl] = 4.00 ; nl++ ;
	  tlist[nl] = RY ; plist[nl] = 90.0 ; nl++ ;
	  tlist[nl] = TX ; plist[nl] = 2.00*(i*2-1) ; nl++ ;
	  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
	  
	  f[onum] = torus; 
	
	  uStart[onum] = 0;		uEnd[onum] = 2*M_PI+0.1,	 uStep[onum] = 0.01;
	  vStart[onum] = 0;		vEnd[onum] = 2*M_PI+0.1,	 vStep[onum] = 0.01;
	
	  onum++ ;


	  // Build larger cylinder ends
	  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
	
	  nl = 0 ;
	  tlist[nl] = SY ; plist[nl] = 10.00 ; nl++ ;
	  tlist[nl] = SZ ; plist[nl] = 10.00 ; nl++ ;
	  tlist[nl] = TX ; plist[nl] = 1.00*(i*2-1) ; nl++ ;
	  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
	  
	  f[onum] = cylinder;
	
	  uStart[onum] = -1;		uEnd[onum] = 1,		 uStep[onum] = 0.01;
	  vStart[onum] = 0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.01;
	
	  onum++ ;



	  // Build doorbays
	  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
	
	  nl = 0 ;
	  tlist[nl] = SY ; plist[nl] = 1.00 ; nl++ ;
	  tlist[nl] = SZ ; plist[nl] = 1.00 ; nl++ ;
	  tlist[nl] = RY ; plist[nl] = 90.0 ; nl++ ;
	  tlist[nl] = TX ; plist[nl] = 2.00*(i*2-1) ; nl++ ;
	  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;

	  
       	  f[onum] = doorbay;
	  
	  uStart[onum] = 0;		uEnd[onum] = 2*M_PI,	 uStep[onum] = 0.01;
	  vStart[onum] = 0;		vEnd[onum] = 1,		 vStep[onum] = 0.01;
	
	  onum++ ;



	  // Build 'spokes'
	  for(int j=0; j<2; ++j) {
		  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
		
		  nl = 0 ;
		  //tlist[nl] = TZ ; plist[nl] = 1.00 ; nl++ ;
		  tlist[nl] = SX ; plist[nl] = 4.00 ; nl++ ;
		  tlist[nl] = SY ; plist[nl] = 4.00 ; nl++ ;
		  tlist[nl] = SZ ; plist[nl] = 4.00 ; nl++ ;
		  tlist[nl] = RY ; plist[nl] = 90.0 ; nl++ ;
		  tlist[nl] = RX ; plist[nl] = 90.00*j ; nl++ ;
		  tlist[nl] = TX ; plist[nl] = 2.00*(i*2-1) ; nl++ ;
		  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
		  
		  f[onum] = cylinder; 
		
		  uStart[onum] = -1;		uEnd[onum] = 1,		 uStep[onum] = 0.01;
		  vStart[onum] = 0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = 0.01;
		
		  onum++ ;

	  }
  }

  // transform all the space-station objects to be in a reasonable position in worldspace
  {
  double tmp_T[4][4], tmp_Ti[4][4];
  nl = 0 ;
  tlist[nl] = RX ; plist[nl] = -55.0 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -30.0 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = 10.0 ; nl++ ;
  M3d_make_movement_sequence_matrix (tmp_T,tmp_Ti,  nl,tlist,plist) ;
  for(int i=0; i<onum; ++i) {
	  M3d_mat_mult(V[i],  tmp_T, V[i]);
	  M3d_mat_mult(Vi[i], Vi[i], tmp_T);
  }
  }











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
  double light_in_world_space[3] = {2,  1, 8};
  
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
		     f,
		     uStart, uEnd, uStep,
		     vStart, vEnd, vStep
		     );






    G_display_image();

    sprintf(filename, "pointcloudimg%04d.xwd", fnum);
    G_save_image_to_file(filename);


    if(G_wait_key() == 'q' || fnum >= max_fnum) {
      break;
    }


    fnum++ ;
  } // end while (1)


  while ('q' != G_wait_key()) {}


}

