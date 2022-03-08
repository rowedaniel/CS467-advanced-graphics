#include "FPToolkit.c"
#include "M3d_matrix_tools.c"

#define M_YON 600

double obmat[100][4][4] ;
double obinv[100][4][4] ;
double color[100][3] ;
int (*grad[100])(double gradient[3], int onum, double intersection[3]);
int    num_objects ;



/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



// define gradient
int circle_grad(double gradient[3], int onum, double intersection[3]) {
  // calculate gradient
  // for circles, gradient is <2x, 2y>
  M3d_mat_mult_pt(gradient, obinv[onum], intersection);
  M3d_vector_mult_const(gradient, gradient, 2);

}


void Draw_ellipsoid (int onum)
{
  int n,i ;
  double t, xyz[3] ;
  double x,y ;

  G_rgb (color[onum][0],color[onum][1],color[onum][2]) ;
  
  n = 1000 ;
  for (i = 0 ; i < n ; i++) {
    t = i*2*M_PI/n ;
    xyz[0] = cos(t) ;
    xyz[1] = sin(t) ;
    xyz[2] = 0 ;
    M3d_mat_mult_pt(xyz, obmat[onum], xyz) ;
    x = xyz[0] ;
    y = xyz[1] ;
    G_point(x,y) ;
  }


}




void Draw_the_scene()
{
  int onum ;
  for (onum = 0 ; onum < num_objects ; onum++) {
    Draw_ellipsoid(onum) ;
  }
}


/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

int solve_quadratic(double a, double b, double c, double t[2])
{
  const double discriminant = b*b-4*a*c;
  if(discriminant< 0) {
	  return 0;
  } else if(discriminant== 0) {
	  t[0] = -b + sqrt(discriminant) / (2 * a);
	  return 1;
  }
  t[0] = (-b + sqrt(discriminant)) / (2 * a);
  t[1] = (-b - sqrt(discriminant)) / (2 * a);
  return 2;
  
}


int get_normal(double normal[3], int onum, double intersection[3]) {
  // get gradient
  double gradient[3];
  grad[onum](gradient, onum, intersection);

  // finally, transform back into object space
  normal[0] = obinv[onum][0][0]*gradient[0] + obinv[onum][1][0]*gradient[1];
  normal[1] = obinv[onum][0][1]*gradient[0] + obinv[onum][1][1]*gradient[1];
  normal[2] = 0;

  M3d_normalize(normal, normal);
}

int ray(double Rsource[3], double Rtip[3], double argb[3])
{
  argb[0] = 1;
  argb[1] = 1;
  argb[2] = 1;

  double t_closest = M_YON + 1;
  int saved_onum = -1;


  // line = start + a*change

  for(int onum = 0; onum < num_objects; ++onum) {
    double tip[3];
    double start[3];
  
    // calculate and start
    M3d_vector_copy(start, Rsource);
    M3d_vector_copy(tip, Rtip);

    // transform line to object space
    M3d_mat_mult_pt(start, obinv[onum], start);
    M3d_mat_mult_pt(tip, obinv[onum], tip);

    // subtract start from change
    double change[3];
    M3d_vector_mult_const(change, start, -1);
    M3d_vector_add(change, change, tip);


    // solve quadratic (only works for objects which started as unit sphere)
    double t[2];
    double a = M3d_dot_product(change, change);
    double b = 2*M3d_dot_product(start, change);
    double c = M3d_dot_product(start, start) -1;
    int num_solutions = solve_quadratic(a, b, c,  t);

    // find closest solution
    for(int i = 0; i < num_solutions; ++i) {
	    if(t[i] < 0) {
		    printf("huh, negative value for solution\n");
		    continue;
	    }
	    if(t[i] < t_closest) {
        // new closest point! Save this one.
	      t_closest = t[i];
        saved_onum = onum;
	      for(int j=0; j<3; ++j) {
	      	argb[j] = color[onum][j];
	      }

	    }

    }

  }

  if(t_closest > M_YON) {
	  return 0;
  }

  // (re) calculate and start
  double change[3];
  M3d_vector_mult_const(change, Rsource, -1);
  M3d_vector_add(change, change, Rtip);

  // finally, add to get back point
  double point[3] = {0,0,0};
  M3d_vector_mult_const(point, change, t_closest);
  M3d_vector_add(point, Rsource, point);
  
  // calculate normal
  double normal[3];
  get_normal(normal, saved_onum, point);




  // finally, draw
  G_rgb(argb[0], argb[1], argb[2]);
  G_circle(Rtip[0], Rtip[1], 2);
  G_line(Rtip[0], Rtip[1], point[0], point[1]);
  G_line(point[0], point[1], point[0] + 20*normal[0], point[1] + 20*normal[1]);
}




int test01()
{
  double vm[4][4], vi[4][4];
  double Tvlist[100];
  int Tn, Ttypelist[100];
  double m[4][4], mi[4][4];
  double Rsource[3];
  double Rtip[3];
  double argb[3] ;

    //////////////////////////////////////////////////////////////////////
    M3d_make_identity(vm) ;    M3d_make_identity(vi) ; // OVERRIDE for 2d
    //////////////////////////////////////////////////////////////////////

    num_objects = 0 ;

    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.0 ;
    color[num_objects][1] = 0.8 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   60   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   25   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  300   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  200   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    grad[num_objects] = circle_grad;
    num_objects++ ; // don't forget to do this

    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  180   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   40   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =   60   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  400   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  550   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    grad[num_objects] = circle_grad;
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.3 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =   75   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   35   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  150   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  360   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  500   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    grad[num_objects] = circle_grad;
    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.5 ;
    color[num_objects][1] = 1.0 ; 
    color[num_objects][2] = 1.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  130   ; Tn++ ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =   30   ; Tn++ ;
    Ttypelist[Tn] = RZ ; Tvlist[Tn] =  -15   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  700   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    grad[num_objects] = circle_grad;
    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.5 ;
    color[num_objects][1] = 0.5 ; 
    color[num_objects][2] = 0.5 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SY ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  100   ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  400   ; Tn++ ;
    Ttypelist[Tn] = TY ; Tvlist[Tn] =  400   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    grad[num_objects] = circle_grad;
    num_objects++ ; // don't forget to do this        
    //////////////////////////////////////////////////////////////

    

    G_rgb(0,0,0) ;
    G_clear() ;

    Draw_the_scene() ;
    
    Rsource[0] =  20 ;  Rsource[1] =  400 ;  Rsource[2] = 0 ;    
    G_rgb(1,0,1) ; G_fill_circle(Rsource[0], Rsource[1], 3) ;
    G_rgb(1,0,1) ; G_line(100,200,  100,600) ;
    
    G_wait_key() ;
    
    double ytip ;
    for (ytip = 200 ; ytip <= 600 ; ytip++) {
      Rtip[0]    = 100 ;  Rtip[1]    = ytip ;  Rtip[2]   = 0  ;    

      G_rgb(1,1,0) ; G_line(Rsource[0],Rsource[1],  Rtip[0],Rtip[1]) ;
      ray (Rsource, Rtip, argb) ; 

      Draw_the_scene() ;
      G_wait_key() ;
    }

    G_rgb(1,1,1) ; G_draw_string("'q' to quit", 50,50) ;
    while (G_wait_key() != 'q') ;
    G_save_image_to_file("2d_Simple_Raytracer.xwd") ;
}




//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////




int main()
{
  G_init_graphics(800,800);
  test01() ;
}
