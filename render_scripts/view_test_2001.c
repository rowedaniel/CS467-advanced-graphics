#include "graph3D.c"



// declare global vars

#define M 100
double numobjects = 0;
int color_type[M];
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
    if(color_type[onum]) {
    	make_graph_inherentrgb_1func(
			uStart[onum], uEnd[onum], uStep[onum],
	       		vStart[onum], vEnd[onum], vStep[onum],
    	       		f[onum], 
               		t,
	       		inherent_rgb[onum]
	       		);
    } else {
    	make_graph_image_1func(
			uStart[onum], uEnd[onum], uStep[onum],
	       		vStart[onum], vEnd[onum], vStep[onum],
    	       		f[onum], 
               		t,
	       		image_IDs[onum], image_size[onum][0], image_size[onum][1]
	       		);
    }

  }
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


int main ()
{
  int onum;

  // stuff for moving each object to the right location
  double V[M][4][4], Vi[M][4][4];
  int nl ;
  int tlist[100] ;
  double plist[100] ;
  double final_trans[4][4], final_transi[4][4];
  double movement_trans[4][4], movement_transi[4][4];
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
  int sphere(double u, double v, double P[3]) {
	  const double a = 1-v*v;
	  if(a<0) { return 0; }
	  const double b = sqrt(a);
  	  P[0] = b*cos(u);
      	  P[1] =  v;
      	  P[2] = b*sin(u);
    	  return 1;
  }

  // "Torus" has center radius 1, outer radius 0.1
  int torus(double u, double v, double P[3]) {
	  const double torus_centerR = 1;
	  const double torus_outerR = 0.1 / 2.0;
	  const double torus_squareness = 1.0 / 3.0;

	  if(v >= 8) { return 0; }


	  double r, z;
	  {
	  const double c = cos(M_PI/2.0 * v);
	  const double s = sin(M_PI/2.0 * v);
	  if(0 <= v && v < 1) {
		  r = 2;
		  z = 2*v-1;
	  } else if (1 <= v && v < 2) {
		  r = 1+s;
		  z = 1-c;
	  } else if(2 <= v && v < 3) {
		  r = 5-2*v;
		  z = 2;
	  } else if(3 <= v && v < 4) {
		  r = -1-c;
		  z = 1-s;
	  } else if(4 <= v && v < 5) {
		  r = -2;
		  z = 9-2*v;
	  } else if(5 <= v && v < 6) {
		  r = -1-s;
		  z = -1+c;
	  } else if(6 <= v && v < 7) {
		  r = -13+2*v;
		  z = -2;
	  } else if(7 <= v && v <= 8) {
		  r = 1+c;
		  z = -1+s;
	  }
	  }


	  P[0] = (torus_centerR + torus_outerR*r) * cos(u);
	  P[1] = (torus_centerR + torus_outerR*r) * sin(u);
	  P[2] = torus_outerR*z;
	  return 1;
  }


  // scaffolding is same as torus, but different parameterization (makes it look cooler)
  int scaffolding(double u, double v, double P[3]) {
	  const double torus_centerR = 1;
	  const double torus_outerR = 0.08;
	  const double torus_squareness = 1.0 / 3.0;
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
  	  P[0] = -u;
  	  P[1] =  rad*cos(v);
  	  P[2] =  rad*sin(v);
	  return 1;
  }

  // The cone has radius of 0.1, centered at origin, 
  // lying on the x-axis from x = 0 to x = 1.
  int cone(double u, double v, double P[3]) {
	  double rad = 0.1;
  	  P[0] = -u;
  	  P[1] =  rad*u*cos(v);
  	  P[2] =  rad*u*sin(v);
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
  double station_highlight[3] = {0.7, 0.3, 0.4};

  double planetres = 0.1;//0.0001;
  double highres = 0.1;//0.001;
  double lowres = 0.2;



  // build planet
  image_IDs[onum] = init_xwd_map_from_file("earth_jeff.xwd");
  if(image_IDs[onum] == -1) { printf("File load failure!\n"); exit(0); }
  error = get_xwd_map_dimensions(image_IDs[onum],image_size[onum]); 
  if(error == -1) { printf("File load failure!\n"); exit(0); }

  nl = 0 ;
  tlist[nl] = RX ; plist[nl] = -45.00   ; nl++ ;
  tlist[nl] = RY ; plist[nl] =  10.00   ; nl++ ;
  tlist[nl] = TX ; plist[nl] =  0.90    ; nl++ ;
  tlist[nl] = TY ; plist[nl] = -0.50    ; nl++ ;
  tlist[nl] = TZ ; plist[nl] = -1.10    ; nl++ ;
  tlist[nl] = SX ; plist[nl] = -140.00  ; nl++ ;
  tlist[nl] = SY ; plist[nl] = -140.00  ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = -140.00  ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  f[onum] = sphere; 

  uStart[onum] = 0;		uEnd[onum] = 2*M_PI, 	 uStep[onum] = planetres;
  vStart[onum] = M_PI/2;	vEnd[onum] = M_PI/2,	 vStep[onum] = planetres;

  onum++ ;


  
  // Build central axis cylinder.
  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
  color_type[onum] = 1;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 2.00 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 6.00 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 6.00 ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  f[onum] = cylinder; 

  uStart[onum] = -1;		uEnd[onum] = 1, 	 uStep[onum] = highres;
  vStart[onum] = 0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = highres;

  onum++ ;


  // Build radial torus
  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
  color_type[onum] = 1;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 4.00 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 4.00 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 4.00 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 90.0 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = 2.00 ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  f[onum] = torus; 

  uStart[onum] = 0;		uEnd[onum] = 2*M_PI;	uStep[onum] = highres;
  vStart[onum] = 0;		vEnd[onum] = 8;		vStep[onum] = highres;
  onum++ ;


  // build scaffolding

  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_highlight[c];}
  color_type[onum] = 1;

  nl = 0 ;
  tlist[nl] = SX ; plist[nl] = 4.00 ; nl++ ;
  tlist[nl] = SY ; plist[nl] = 4.00 ; nl++ ;
  tlist[nl] = SZ ; plist[nl] = 4.00 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = 90.0 ; nl++ ;
  tlist[nl] = TX ; plist[nl] = -2.00 ; nl++ ;
  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
  
  f[onum] = scaffolding; 

  uStart[onum] = 0;		uEnd[onum] = 2*M_PI;	uStep[onum] = lowres;
  vStart[onum] = 0;		vEnd[onum] = 8;		vStep[onum] = lowres;

  onum++ ;

  for(int i=0; i<4; ++i) {
	  // built sections of scaffold
	  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
  	  color_type[onum] = 1;

	  nl = 0 ;
	  tlist[nl] = SX ; plist[nl] = 4.00 ; nl++ ;
	  tlist[nl] = SY ; plist[nl] = 4.00 ; nl++ ;
	  tlist[nl] = SZ ; plist[nl] = 4.00 ; nl++ ;
	  tlist[nl] = RY ; plist[nl] = 90.0 ; nl++ ;
	  tlist[nl] = RX ; plist[nl] = i*90.0 ; nl++ ;
	  tlist[nl] = TX ; plist[nl] = -2.00 ; nl++ ;
	  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
	  
	  f[onum] = torus; 

	  uStart[onum] = -M_PI/10.0;	uEnd[onum] = M_PI/10.0;	uStep[onum] = highres;
	  vStart[onum] = 0;		vEnd[onum] = 8;		vStep[onum] = highres;

	  onum++ ;
  }



  for(int i=0; i<2; ++i) {


	  // Build larger cylinder ends
	  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
  	  color_type[onum] = 1;
	
	  nl = 0 ;
	  tlist[nl] = SY ; plist[nl] = 10.00 ; nl++ ;
	  tlist[nl] = SZ ; plist[nl] = 10.00 ; nl++ ;
	  tlist[nl] = SX ; plist[nl] = 0.50 ; nl++ ;
	  tlist[nl] = TX ; plist[nl] = 2.00*(i*2-1) ; nl++ ;
	  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
	  
	  f[onum] = cylinder;
	
	  uStart[onum] = -1;		uEnd[onum] = 1,		 uStep[onum] = highres;
	  vStart[onum] = 0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = highres;
	
	  onum++ ;

	  // build cone connecters for larger ends
	  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
  	  color_type[onum] = 1;
	
	  nl = 0 ;
	  tlist[nl] = SY ; plist[nl] = 10.00 ; nl++ ;
	  tlist[nl] = SZ ; plist[nl] = 10.00 ; nl++ ;
	  tlist[nl] = SX ; plist[nl] = 1.00  ; nl++ ;
	  tlist[nl] = TX ; plist[nl] = 0.50  ; nl++ ;
	  tlist[nl] = RY ; plist[nl] = 180.0*i ; nl++ ;
	  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
	  
	  f[onum] = cone;
	
	  uStart[onum] = -1;		uEnd[onum] = 0,		 uStep[onum] = highres;
	  vStart[onum] = 0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = highres;
	
	  onum++ ;



	  // Build doorbays
	  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
  	  color_type[onum] = 1;
	
	  nl = 0 ;
	  tlist[nl] = SY ; plist[nl] = 1.00 ; nl++ ;
	  tlist[nl] = SX ; plist[nl] = 1.00 ; nl++ ;
	  tlist[nl] = RY ; plist[nl] = -90.0; nl++ ;
	  tlist[nl] = TX ; plist[nl] = 2.50 ; nl++ ;
	  tlist[nl] = RY ; plist[nl] = 180.0*i; nl++ ;
	  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;

	  
       	  f[onum] = doorbay;
	  
	  uStart[onum] = 0;		uEnd[onum] = 2*M_PI,	 uStep[onum] = highres/10;
	  vStart[onum] = 0;		vEnd[onum] = 1,		 vStep[onum] = highres/10;
	
	  onum++ ;



	  // Build 'spokes'
	  for(int j=0; j<4; ++j) {
		  for(int c=0;c<3;++c) {inherent_rgb[onum][c]=station_color[c];}
  		  color_type[onum] = 1;
		
		  nl = 0 ;
		  //tlist[nl] = TZ ; plist[nl] = 1.00 ; nl++ ;
		  tlist[nl] = SX ; plist[nl] = 1.50 ; nl++ ;
		  tlist[nl] = SY ; plist[nl] = 4.00 ; nl++ ;
		  tlist[nl] = SZ ; plist[nl] = 4.00 ; nl++ ;
		  tlist[nl] = RY ; plist[nl] = 90.0 ; nl++ ;
		  tlist[nl] = TZ ; plist[nl] = 2.50 ; nl++ ;
		  tlist[nl] = RX ; plist[nl] = 90.00*j ; nl++ ;
		  tlist[nl] = TX ; plist[nl] = 2.00*(i*2-1) ; nl++ ;
		  M3d_make_movement_sequence_matrix (V[onum],Vi[onum],  nl,tlist,plist) ;
		  
		  f[onum] = cylinder; 
		
		  uStart[onum] = -1;		uEnd[onum] = 1,		 uStep[onum] = highres;
		  vStart[onum] = 0;		vEnd[onum] = 2*M_PI,	 vStep[onum] = highres;
		
		  onum++ ;

	  }
  }

  // transform all the space-station objects to be in a reasonable position in worldspace
  nl = 0 ;
  tlist[nl] = RX ; plist[nl] = -55.0 ; nl++ ;
  tlist[nl] = RY ; plist[nl] = -50.0 ; nl++ ;
  tlist[nl] = RZ ; plist[nl] = 10.0 ; nl++ ;
  M3d_make_movement_sequence_matrix (final_trans,final_transi,  nl,tlist,plist) ;












  //---------------------------------------------------------------------
  numobjects = onum ;
  //---------------------------------------------------------------------


  // object movement stuff
  double final_V[M][4][4];

  // eye stuff
  double eye[3], coi[3], up[3] ;
  double v[4][4], vi[4][4];
  int i;

  int fnum ;
  int max_fnum;
  double t ;

  fnum = 0 ;
  max_fnum = 360;


  // light model setup
  double light_in_world_space[3] = {5,  1, 0};
  
  AMBIENT = 0.2 ;
  MAX_DIFFUSE = 0.5 ;
  SPECPOW = 30 ;

  // file saving setup
  char filename[100];


  while (1) {

    t = 0.003333*fnum ;

    eye[0] = 15-t*10; //15(2*M_PI*t) ; 
    eye[1] =  0; //6*t ; 
    eye[2] =  5; //7*sin(2*M_PI*t) ; 

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
    
    
    // update space station rotation
    nl = 0 ;
    tlist[nl] = RX ; plist[nl] = t*100 ; nl++ ;
    M3d_make_movement_sequence_matrix (movement_trans,movement_transi,  nl,tlist,plist) ;

    // copy transformation over for planet (it doesn't move)
    M3d_copy_mat(final_V[0], V[0]);
    // update transformation for moving objects (not obj 0--the planet)
    for(int i=1; i<onum; ++i) {
	    M3d_mat_mult(final_V[i], final_trans, movement_trans);
	    M3d_mat_mult(final_V[i], final_V[i], V[i]);
    }

    make_z_buff();
    draw_all_objects(v, final_V,
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

