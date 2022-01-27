
#include "FPToolkit.c"
#include "M3d_matrix_tools.c"
#include "time.h"

// stickfigure initially designed as centered for a 400x400 window :
double x[13] = {175,225,225,300,225,225,250,200,150,175,175,100,175} ;
double y[13] = {300,300,250,225,225,200,100,175,100,200,225,225,250} ;
double z[13] = {0,0,0,0,0,0,0,0,0,0,0,0,0} ;
       // z[] values unimportant but should NOT be left uninitialized
       // as nan values WILL propagate through
int n = 13 ;



int main(int argc, char **argv) 
{

 if (argc != 3) {
    printf("usage : pgm_name   window_size  microseconds(30000)\n") ;
    exit(0) ;
 }
  
 double winsize = atof(argv[1]) ;
 double u = atof(argv[2]) ; 

 G_init_graphics(winsize,winsize) ;


 // the original design was for a 400x400
 // window and the object is centered on 200,200
 // so we recenter it and make it larger
 // (you get to do this ... use the
 // M3d_make_movement_sequence_matrix  function :)

 // .....


 const double rotspeed = 2.0;
 const double scalespeed = 0.95;

 double v[4][4], vi[4][4];
 int i ;
 int mtype[100] ;
 double mparam[100] ;

 // rescale stick-figure to match the screensize
 i = 0 ;
 mtype[i] = TX ;  mparam[i] =  -200      ; i++ ;
 mtype[i] = TY ;  mparam[i] =  -200      ; i++ ;
 mtype[i] = SX ;  mparam[i] =  winsize / 400.0 ; i++ ;
 mtype[i] = SY ;  mparam[i] =  winsize / 400.0 ; i++ ;
 mtype[i] = TX ;  mparam[i] =  winsize / 2.0   ; i++ ;
 mtype[i] = TY ;  mparam[i] =  winsize / 2.0   ; i++ ;

 M3d_make_movement_sequence_matrix(v,vi,  i,mtype,mparam) ;
 M3d_mat_mult_points(x, y, z, v, x, y, z, n);


 // now make the movie the rotates and shrinks about the center :
 i = 0 ;
 mtype[i] = TX ;  mparam[i] =  -winsize/2.0 ; i++ ;
 mtype[i] = TY ;  mparam[i] =  -winsize/2.0 ; i++ ;
 mtype[i] = RZ ;  mparam[i] =  rotspeed     ; i++ ;
 mtype[i] = SX ;  mparam[i] =  scalespeed   ; i++ ;
 mtype[i] = SY ;  mparam[i] =  scalespeed   ; i++ ;
 mtype[i] = TX ;  mparam[i] =  winsize/2.0  ; i++ ;
 mtype[i] = TY ;  mparam[i] =  winsize/2.0  ; i++ ;

 M3d_make_movement_sequence_matrix(v,vi,  i,mtype,mparam) ;

 printf("2nd mat:\n");
 M3d_print_mat(v);

 // .....
 
 G_rgb(0,0,0);
 G_clear();
 G_rgb(1,1,1);
 G_fill_polygon(x,y,n);
 G_wait_key();
 

 printf("frametime: %lf\n", u);
 for(i=0; i < (360/rotspeed); ++i) { 

   G_rgb(0,0,0);
   G_clear();

   G_rgb(0.5, 0.4, 0.2);
   G_fill_polygon(x, y, n);
   M3d_mat_mult_points(x, y, z, v, x, y, z, n);

   G_display_image();

   /*
   if(G_no_wait_key() == 'q') {
     break;
   }
   */

   usleep(u);

 }

 G_rgb(1, 0, 0);
 G_fill_circle(winsize/2.0, winsize/2.0, winsize/4.0);
 G_wait_key();

}

