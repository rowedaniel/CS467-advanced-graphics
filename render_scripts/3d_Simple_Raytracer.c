#include "FPToolkit.c"
#include "M3d_matrix_tools.c"

#define M_YON 10000

double obmat[100][4][4] ;
double obinv[100][4][4] ;
double color[100][3] ;
int (*grad[100])(double gradient[3], int onum, double intersection[3]);
void (*draw[100])(int onum);
double (*intersection[100])(double start[3], double change[3]);
int    num_objects ;



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



// ================ Sphere stuff ===============
int sphere_grad(double gradient[3], int onum, double intersection[3]) {
  // for spheres, gradient is <2x, 2y, 2z>
  M3d_mat_mult_pt(gradient, obinv[onum], intersection);
  M3d_vector_mult_const(gradient, gradient, 2);
}

double sphere_intersection(double start[3], double change[3])
{
  // solve quadratic (only works for objects which started as unit sphere)
  double a = M3d_dot_product(change, change);
  double b = 2*M3d_dot_product(start, change);
  double c = M3d_dot_product(start, start) -1;

  double t[2];
  double t_closest = M_YON + 1;
  int num_solutions = solve_quadratic(a, b, c,  t);

  // find closest solution
  for(int i = 0; i < num_solutions; ++i) {
    if(t[i] < 0) {
      continue;
    }
    if(t[i] < t_closest) {
      t_closest = t[i];
    }
  }

  return t_closest;
}

// ==========================================================================


// ================ plane stuff ===============
int plane_grad(double gradient[3], int onum, double intersection[3]) {
  gradient[0] = 0;
  gradient[1] = 0;
  gradient[2] = 1;
}
double plane_intersection(double start[3], double change[3]) {
  if(change[2] == 0) {
    return M_YON + 1;
  }
  double t = -start[2] / change[2];

  double x = start[0] + t * change[0];
  double y = start[1] + t * change[1];
  if(-1 > x || x > 1 ||
     -1 > y || y > 1) {
    return M_YON + 1;
  }
  return t;
}




/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////



int get_normal(double normal[3], int onum, double intersection[3]) {
  // get gradient
  double gradient[3];
  grad[onum](gradient, onum, intersection);

  // finally, transform back into object space
  double inv_T[4][4];
  M3d_transpose(inv_T, obinv[onum]);
  M3d_mat_mult_pt(normal, inv_T, gradient);

  M3d_normalize(normal, normal);
}

int ray_recursive(double Rsource[3], double Rtip[3], double argb[3], int n)
{

  double t_closest = M_YON + 1;
  int saved_onum = -1;


  // line = start + t*change
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

    double t = intersection[onum](start, change);
    if(t < t_closest) {
      // new closest point! Save this one.
      t_closest = t;
      saved_onum = onum;
    }


  }

  if(t_closest > M_YON) {
	  return 0;
  }


  // set color to this object's color:
  // TODO: implement lightmodel here
  // TODO: implement reflection here
  if(n == 0) {
    for(int j=0; j<3; ++j) {
      argb[j] = color[saved_onum][j];
    }
  }



  // recurse!
  if (n > 0) {
    // get a point on the opposite side of Rtip as Rsource
    // first, get the change = Rtip - Rsource
    double change[3];
    M3d_vector_mult_const(change, Rsource, -1);
    M3d_vector_add(change, change, Rtip);

    // then get point = Rsource + t(change)
    double point[3] = {0,0,0};
    // reduce t_closest slightly to avoid z-fighting
    t_closest *= 0.999;
    M3d_vector_mult_const(point, change, t_closest);
    M3d_vector_add(point, Rsource, point);

    // calculate normal
    double normal[3];
    get_normal(normal, saved_onum, point);

    // calculate reflection angle
    // get transformation matrix to reflect across plane defined by normal vector
    double look[3], reflection[3];
    // look = intersection - source
    M3d_vector_mult_const(look, Rsource, -1);
    M3d_vector_add(look, look, point);
    M3d_normalize(look, look);

    // reflection = look - 2(look * normal)normal
    M3d_vector_mult_const(reflection, normal, -2*M3d_dot_product(look, normal));
    M3d_vector_add(reflection, reflection, look);

    if(saved_onum == 0) {
      printf("old change vector:\n");
      M3d_print_vector(change);
      printf("new look vector:\n");
      M3d_print_vector(look);
    }

    // new tip = reflection + intersection
    double new_tip[3];
    M3d_vector_add(new_tip, reflection, point);


    ray_recursive(point, new_tip, argb, n-1);
  }



  


}

int ray(double Rsource[3], double Rtip[3], double argb[3])
{
  ray_recursive(Rsource, Rtip, argb, 1);
}


int screen_to_coords(double point[3], double screen_pos[2]) {
  point[0] = screen_pos[0] * 2 / 800 - 1;
  point[1] = screen_pos[1] * 2 / 800 - 1;
  point[2] = 1;
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
    color[num_objects][0] = 1.0 ;
    color[num_objects][1] = 0.3 ; 
    color[num_objects][2] = 0.0 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = RY ; Tvlist[Tn] =  45  ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  4   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    grad[num_objects] = plane_grad;
    intersection[num_objects] = plane_intersection;
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.2 ;
    color[num_objects][1] = 0.6 ; 
    color[num_objects][2] = 0.3 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = RY ; Tvlist[Tn] =  90  ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  -2   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  4   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    grad[num_objects] = plane_grad;
    intersection[num_objects] = plane_intersection;
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    color[num_objects][0] = 0.3 ;
    color[num_objects][1] = 0.2 ; 
    color[num_objects][2] = 0.6 ;
	
    Tn = 0 ;
    Ttypelist[Tn] = SX ; Tvlist[Tn] =  0.1  ; Tn++ ;
    Ttypelist[Tn] = RY ; Tvlist[Tn] =  90  ; Tn++ ;
    Ttypelist[Tn] = TX ; Tvlist[Tn] =  3   ; Tn++ ;
    Ttypelist[Tn] = TZ ; Tvlist[Tn] =  4   ; Tn++ ;
	
    M3d_make_movement_sequence_matrix(m, mi, Tn, Ttypelist, Tvlist);
    M3d_mat_mult(obmat[num_objects], vm, m) ;
    M3d_mat_mult(obinv[num_objects], mi, vi) ;

    grad[num_objects] = plane_grad;
    intersection[num_objects] = plane_intersection;
    num_objects++ ; // don't forget to do this
    //////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////

    

    Rsource[0] =  0 ;  Rsource[1] =  0 ;  Rsource[2] = 0 ;    

    G_rgb(0,0,0) ;
    G_clear() ;

    // go through a first random pass to draw some points
    for(int i=0; i<800; ++i) {
      for(int j=0; j<800; ++j) {
        double screen_pos[2];
        screen_pos[0] = i * 1.0;
        screen_pos[1] = j * 1.0;
        screen_to_coords(Rtip, screen_pos);

        argb[0] = 1;
        argb[1] = 1;
        argb[2] = 1;
        ray (Rsource, Rtip, argb) ; 
        G_rgb(argb[0], argb[1], argb[2]);
        G_point(screen_pos[0], screen_pos[1]);
      }
    }

    while(1)
    {



      double p[2];
      G_wait_click(p);
      if(p[1] < 50) { break; }

      screen_to_coords(Rtip, p);



      // redraw screen
      //G_rgb(0,0,0) ;
      //G_clear() ;
  
      G_rgb(1,1,1) ; G_draw_string("click here to quit", 50,50) ;


      // draw point
      argb[0] = 1;
      argb[1] = 1;
      argb[2] = 1;
      ray (Rsource, Rtip, argb) ; 
      G_rgb(argb[0], argb[1], argb[2]);
      G_point(p[0], p[1]);

    }

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
